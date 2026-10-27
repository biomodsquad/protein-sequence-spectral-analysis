"""Build a self-contained Plotly viewer from the completed atlas (optional dependency)."""

import argparse
import html
import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .analyze import COLORS


def main():
    """Export hoverable charts, identifier lookup, aliases, and whole-geometry neighbors."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    args = parser.parse_args()
    root = args.results
    vectors = pd.read_csv(root/'unique_vector_layouts.csv')
    whole = pd.read_csv(root/'whole_representations.csv')
    catalogue = pd.read_csv(root/'vector_catalogue.csv')
    cka = pd.read_csv(root/'whole_cka_equal_channels.csv', index_col=0)
    aliases = catalogue.groupby('unique_id').id.apply(list).to_dict()
    vectors['alias_names'] = vectors.unique_id.map(lambda key: ', '.join(aliases[key]))
    vectors['hover'] = vectors.apply(lambda r: html.escape(
        f'{r.id} | cluster {r.cluster_r095} | {r.aliases} aliases: {r.alias_names[:450]}'), axis=1)
    whole['hover'] = whole.apply(lambda r: html.escape(
        f'{r.id}: {r.description} | effective rank {r.effective_rank:.2f}'), axis=1)
    fig = make_subplots(rows=2, cols=2, subplot_titles=(
        'Vectors: sign-invariant UMAP', 'Vectors: projector PCA',
        'Whole representations: CKA-space UMAP', 'Whole representations: CKA-space PCA'))
    layouts = [(vectors, 'umap_n15_seed0_1', 'umap_n15_seed0_2', 1, 1),
               (vectors, 'projector_PC1', 'projector_PC2', 1, 2),
               (whole, 'UMAP1', 'UMAP2', 2, 1), (whole, 'PC1', 'PC2', 2, 2)]
    for panel, (table, x, y, row, col) in enumerate(layouts):
        for family, color in COLORS.items():
            rows = table[table.family == family]
            reference = family in ('BLOSUM100', 'Kidera')
            fig.add_trace(go.Scattergl(x=rows[x], y=rows[y], mode='markers', name=family,
                                      legendgroup=family, showlegend=panel == 0,
                                      text=rows.hover, hovertemplate='%{text}<extra></extra>',
                                      marker=dict(color=color, size=10 if reference else 5,
                                                  symbol='star' if reference else 'circle',
                                                  opacity=1 if reference else .6)), row=row, col=col)
    highlights = []
    for row, col in [(1, 1), (1, 2), (2, 1), (2, 2)]:
        highlights.append(len(fig.data))
        fig.add_trace(go.Scatter(x=[], y=[], mode='markers', showlegend=False,
                                marker=dict(color='red', size=19, symbol='circle-open', line=dict(width=3))),
                      row=row, col=col)
    fig.update_layout(height=1000, template='plotly_white', margin=dict(l=45,r=25,t=55,b=25))
    plot = fig.to_html(full_html=False, include_plotlyjs=True, div_id='atlas')
    lookup = {}
    vector_rows = vectors.set_index('unique_id')
    for r in catalogue.itertuples():
        point = vector_rows.loc[r.unique_id]
        lookup[r.id] = {'u':[float(point.umap_n15_seed0_1),float(point.umap_n15_seed0_2)],
                        'p':[float(point.projector_PC1),float(point.projector_PC2)],
                        'text':f'{r.id}: {r.description}; aliases: '+', '.join(aliases[r.unique_id])}
    whole_lookup = {}
    for r in whole.itertuples():
        scores = cka.loc[r.id].drop(r.id).nlargest(8)
        whole_lookup[r.id] = {'u':[r.UMAP1,r.UMAP2],'p':[r.PC1,r.PC2],
                              'text':r.description, 'neighbors':[[k,float(v)] for k,v in scores.items()]}
    options = ''.join(f'<option value="{html.escape(key,quote=True)}">' for key in lookup)
    choices = ''.join(f'<option>{html.escape(key)}</option>' for key in whole_lookup)
    header = '''<!doctype html><html><head><meta charset="utf-8"><title>Protspect embedding atlas</title>
<style>body{font:16px system-ui;max-width:1500px;margin:auto;padding:20px;color:#203040}
input,select,button{font:inherit;padding:8px}input{width:430px}#detail{white-space:pre-wrap;
background:#f0f4f7;padding:15px;border-radius:8px;overflow-wrap:anywhere}p{max-width:1100px}</style>
</head><body><h1>Protspect embedding atlas</h1><p>6,051 valid channel instances → 3,826 unique centered
shapes; 166 whole representations. Hover to inspect, drag to zoom, and click legend entries to hide/show
families. Family colors choose one source for shared shapes; all aliases remain available below.</p>
<p>Vector distances use 1−r² after centering/unit scaling (sign-invariant). Whole similarities use
linear CKA after channel centering/equalization. Plotted distances are approximations, not clustering
criteria or predictive-performance measurements. PCA axes retain only part of the full variation.</p>'''
    header = (header.replace('6,051', f'{len(catalogue):,}')
              .replace('3,826', f'{len(vectors):,}').replace('166 whole', f'{len(whole)} whole'))
    controls = f'''<label>Locate a channel: <input id="query" list="ids" placeholder="e.g. Kidera:KF4 or KIDA850101"></label>
<datalist id="ids">{options}</datalist><button id="find">Locate</button><br><br>
<label>Whole representation: <select id="whole">{choices}</select></label><button id="compare">Compare</button>
<pre id="detail">Choose an identifier to inspect aliases or nearest whole representations.</pre>'''
    script = '''<script>
const lookup=LOOKUP, whole=WHOLE, selected=HIGHLIGHTS;
function mark(points,start){points.forEach((p,i)=>Plotly.restyle('atlas',{x:[[p[0]]],y:[[p[1]]]},[selected[start+i]]));}
document.getElementById('find').onclick=()=>{const key=document.getElementById('query').value;
 const p=lookup[key]; if(!p){document.getElementById('detail').textContent='Identifier not found.';return;}
 mark([p.u,p.p],0);document.getElementById('detail').textContent=p.text;};
document.getElementById('compare').onclick=()=>{const key=document.getElementById('whole').value,p=whole[key];
 mark([p.u,p.p],2);document.getElementById('detail').textContent=key+': '+p.text+'\\n\\nNearest whole geometries (equal-channel CKA):\\n'+p.neighbors.map(r=>r[0]+'  '+r[1].toFixed(4)).join('\\n');};
</script></body></html>'''.replace('LOOKUP',json.dumps(lookup).replace('</','<\\/')).replace(
        'WHOLE',json.dumps(whole_lookup).replace('</','<\\/')).replace('HIGHLIGHTS',json.dumps(highlights))
    path = root/'interactive_atlas.html'
    path.write_text(header+controls+plot+script,encoding='utf-8')
    print(path, 'bytes', path.stat().st_size)
    # Validate source/alias/coordinate coverage independent of rendering.
    assert len(lookup) == len(catalogue) and len(whole_lookup) == len(whole)
    assert all(np.isfinite(item['u']+item['p']).all() for item in lookup.values())


if __name__ == '__main__':
    main()
