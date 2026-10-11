#!/usr/bin/env python3
"""Reproducible BPQ-style molecular case PDFs (no invented molecules or scores).

Inputs MUST contain actual source + generated SMILES, and benchmark-calculated
property scores. Missing structures are intentionally not visualized.
"""
import json, sys
from pathlib import Path
from io import BytesIO
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor

INK=HexColor("#162839")
BLUE=HexColor("#287FAE")
GREEN=HexColor("#07975B")
GREY=HexColor("#617389")
RULE=HexColor("#D7E1E8")
W, PANEL, GAP=755, 197, 22

def txt(c,x,y,value,size=10,bold=False,color=INK,align="left"):
    c.setFillColor(color)
    c.setFont("Helvetica-Bold" if bold else "Helvetica",size)
    value=str(value)
    if align=="right": c.drawRightString(x,y,value)
    elif align=="center": c.drawCentredString(x,y,value)
    else: c.drawString(x,y,value)

def molecule(c,smiles,x,y,w,h):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError("Invalid SMILES: "+smiles)
    rdDepictor.Compute2DCoords(mol)
    drawer=rdMolDraw2D.MolDraw2DSVG(480,220)
    drawer.drawOptions().bondLineWidth=2
    drawer.drawOptions().padding=.075
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    graph=svg2rlg(BytesIO(drawer.GetDrawingText().encode("utf-8")))
    scale=min(w/graph.width,h/graph.height)
    c.saveState()
    c.translate(x+(w-graph.width*scale)/2,y+(h-graph.height*scale)/2)
    c.scale(scale,scale)
    renderPDF.draw(graph,c,0,0)
    c.restoreState()

def panel(c,case,out,ox,oy):
    c.setStrokeColor(RULE)
    c.setLineWidth(.8)
    c.roundRect(ox+.5,oy+.5,W-1,PANEL-1,7,stroke=1,fill=0)
    txt(c,ox+12,oy+PANEL-19,"SOURCE MOLECULE",9,True,GREY)
    txt(c,ox+W/2,oy+PANEL-19,f"SOR {out['SOR']} / SSOR {out['SSOR']}",9,True,GREEN,"center")
    txt(c,ox+W-12,oy+PANEL-19,out["name"],10,True,BLUE,"right")
    molecule(c,case["source"]["smiles"],ox+12,oy+37,315,129)
    molecule(c,out["smiles"],ox+W-327,oy+37,315,129)
    cx=ox+W/2;cy=oy+102
    c.setStrokeColor(GREY);c.setLineWidth(1.4)
    p=c.beginPath();p.moveTo(cx-24,cy-5);p.lineTo(cx+12,cy-5)
    p.lineTo(cx+12,cy-13);p.lineTo(cx+28,cy);p.lineTo(cx+12,cy+13)
    p.lineTo(cx+12,cy+5);p.lineTo(cx-24,cy+5);p.close()
    c.drawPath(p,stroke=1,fill=0)
    src=case["source"]["metrics"];m=out["metrics"]
    txt(c,ox+12,oy+23,"Source: "+"   ".join(f"{k} {v:.4f}" for k,v in src.items()),9,color=GREY)
    diffs=[f"{k} {m[k]:.4f} ({m[k]-src[k]:+.4f})" for k in ["AMP","BBBP"] if k in src and k in m]
    txt(c,ox+12,oy+9,"Generated: "+"   ".join(diffs),10,True,GREEN)
    rest="   ".join(f"{k} {m[k]:.4f}" for k in ["MUTAG","Sim"] if k in m)
    if rest: txt(c,ox+W-12,oy+9,rest,9,color=GREY,align="right")

def render(case,out_dir):
    if not case.get("source",{}).get("smiles"): raise ValueError("Missing source SMILES")
    outputs=case.get("candidates",[])
    if not outputs or any(not x.get("smiles") for x in outputs):
        raise ValueError("Missing actual model-generated SMILES for "+case["id"])
    height=len(outputs)*PANEL+(len(outputs)-1)*GAP+42
    output=out_dir/("cmoral_case_"+case["id"]+".pdf")
    c=canvas.Canvas(str(output),pagesize=(W,height),pageCompression=1)
    c.setTitle(case["title"]+" - source and model generated molecular structures")
    txt(c,12,height-17,case["title"],13,True)
    txt(c,W-12,height-17,case["backbone"]+" | greedy candidate",9,color=GREY,align="right")
    for j,out in enumerate(outputs):
        panel(c,case,out,0,height-42-(j+1)*PANEL-j*GAP)
    c.showPage();c.save()
    print("Generated",output,output.stat().st_size,"bytes")

if __name__=="__main__":
    root=Path(__file__).resolve().parents[1]
    data=json.loads((root/"latex/case_molecule_data.json").read_text(encoding="utf-8"))
    figs=root/"figs"
    for case in data["cases"]: render(case,figs)
    pending=data.get("pending_case_structures",[])
    if pending: print("Not rendered (original generated SMILES missing):",", ".join(pending),file=sys.stderr)
