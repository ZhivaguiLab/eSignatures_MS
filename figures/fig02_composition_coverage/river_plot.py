import csv, collections
import plotly.graph_objects as go

rows = list(csv.DictReader(open('large_data_ames_iarc_ChemGroup_v3.txt'), delimiter='\t'))

def norm_type(t):
    t = t.lower().replace('-', '')
    return 'in vivo' if t == 'invivo' else 'in vitro'

def norm_seq(s):
    return 'WGS' if s.strip().lower() == 'genome' else 'WES'

def norm_species(s):
    s = s.strip()
    return s if s in ('C. elegans', 'Human', 'Mouse') else 'Others'

def freq(r):
    try: return float(r['Freq'])
    except: return 0.0

# weighted records: (species, type, sequencing, weight)
recs = [(norm_species(r['species']), norm_type(r['type']), norm_seq(r['sequencing']), freq(r))
        for r in rows]
TOTAL = sum(w for *_, w in recs)

specs = ['C. elegans', 'Mouse', 'Human', 'Others']
types = ['in vivo', 'in vitro']
seqs  = ['WGS', 'WES']

# v2 palette: C. elegans gold, Human green, Mouse unchanged
spec_colors = {'C. elegans': '#D4A23A', 'Mouse': '#C97B84', 'Human': '#5C9A6B', 'Others': '#9AA0A6'}
type_colors = {'in vivo': '#3D5A80', 'in vitro': '#98C1D9'}
seq_colors  = {'WGS': '#C1666B', 'WES': '#E9C46A'}
COLOR = {**spec_colors, **type_colors, **seq_colors}

# column x positions (close together -> short flows)
COLX = {'spec': 0.10, 'type': 0.50, 'seq': 0.90}
PAD = 0.025  # vertical gap fraction between nodes in a column

def stack(items, vals):
    n = len(items)
    span = 1 - PAD * (n - 1)
    pos = {}
    y = 0.0
    for it in items:
        h = vals[it] / TOTAL * span
        pos[it] = (y, y + h, y + h / 2)  # top, bottom, center (0 = top)
        y += h + PAD
    return pos

def build(mode):  # mode: 'count' or 'pct'
    tot_sp = {sp: sum(w for s, _, _, w in recs if s == sp) for sp in specs}
    tot_ty = {ty: sum(w for _, t, _, w in recs if t == ty) for ty in types}
    tot_sq = {sq: sum(w for _, _, q, w in recs if q == sq) for sq in seqs}

    ps = stack(specs, tot_sp)
    pt = stack(types, tot_ty)
    pq = stack(seqs, tot_sq)

    order = specs + types + seqs
    idx = {n: i for i, n in enumerate(order)}
    nx = [COLX['spec']] * len(specs) + [COLX['type']] * len(types) + [COLX['seq']] * len(seqs)
    ny = ([ps[s][2] for s in specs] + [pt[t][2] for t in types] + [pq[q][2] for q in seqs])
    colors = [COLOR[n] for n in order]

    def hx(name, a=0.42):
        c = COLOR[name].lstrip('#'); r, g, b = (int(c[i:i+2], 16) for i in (0, 2, 4))
        return f'rgba({r},{g},{b},{a})'

    src, tgt, val, lcol = [], [], [], []
    c1 = collections.Counter(); c2 = collections.Counter()
    for sp, ty, sq, w in recs:
        c1[(sp, ty)] += w; c2[(ty, sq)] += w
    for (sp, ty), w in c1.items():
        src.append(idx[sp]); tgt.append(idx[ty]); val.append(w); lcol.append(hx(sp))
    for (ty, sq), w in c2.items():
        src.append(idx[ty]); tgt.append(idx[sq]); val.append(w); lcol.append(hx(ty))

    fig = go.Figure(go.Sankey(
        arrangement='fixed',
        node=dict(label=[''] * len(order), x=nx, y=ny, pad=20, thickness=18,
                  color=colors, line=dict(color='rgba(255,255,255,0)', width=0)),
        link=dict(source=src, target=tgt, value=val, color=lcol,
                  line=dict(width=0))))

    # ---- manual annotations (no white background) ----
    def num(name, tot):
        if mode == 'pct':
            return f'{tot/TOTAL*100:.1f}%'
        return f'{int(round(tot))}'

    ann = []
    def add(name, col_x, cy, side, tot):
        # side: 'left' label sits left of node, 'right' right of node, 'mid' above node
        if side == 'left':
            x, xa = col_x - 0.035, 'right'
        elif side == 'right':
            x, xa = col_x + 0.035, 'left'
        ann.append(dict(x=x, y=1 - cy, xref='paper', yref='paper', xanchor=xa,
                        yanchor='middle',
                        align=('right' if side == 'left' else 'left'),
                        showarrow=False,
                        text=f'<b>{name}</b><br><span style="font-size:11px;color:#666">{num(name, tot)}</span>',
                        font=dict(size=13, color='#2b2b2b',
                                  family='Helvetica Neue, Arial')))

    for s in specs:
        add(s, COLX['spec'], ps[s][2], 'left', tot_sp[s])
    for q in seqs:
        add(q, COLX['seq'], pq[q][2], 'right', tot_sq[q])
    # middle column labels centered vertically on each node's bar
    for t in types:
        ann.append(dict(x=COLX['type'], y=1 - pt[t][2], xref='paper', yref='paper',
                        xanchor='center', yanchor='middle', showarrow=False, align='center',
                        text=f'<b>{t}</b><br><span style="font-size:11px;color:#666">{num(t, tot_ty[t])}</span>',
                        font=dict(size=13, color='#2b2b2b', family='Helvetica Neue, Arial')))

    fig.update_layout(
        title=dict(text='Sample Composition', x=0.015, xanchor='left', y=0.97,
                   font=dict(size=20, color='#222', family='Helvetica Neue, Arial')),
        annotations=ann,
        paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
        width=620, height=540, margin=dict(l=80, r=80, t=75, b=30))
    return fig

for mode, name in [('count', 'counts'), ('pct', 'pct')]:
    fig = build(mode)
    fig.write_image(f'river_plot_{name}_v2.pdf', scale=2)
    fig.write_html(f'river_plot_{name}_v2.html')
    print('wrote river_plot_%s_v2.pdf' % name)
