"""Genera datos, figuras, tablas LaTeX y PDF desde una única implementación."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

os.environ.setdefault("SOURCE_DATE_EPOCH", "1791331200")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from .modelo import (A, B, C, C_RK, R, T, X_MAX, aluminio, hierro, energia,
                     diferencia, coeficientes_gamma, limite_uno_por_ciento)

RAIZ = Path(__file__).resolve().parents[2]
SALIDA = RAIZ / "output"


def escribir_csv(nombre, cabecera, filas):
    with (SALIDA / "datos" / nombre).open("w", encoding="utf-8", newline="") as f:
        escritor = csv.writer(f)
        escritor.writerow(cabecera)
        escritor.writerows(filas)


def tabla_tex(nombre, cabecera, filas, formato):
    texto = "\\begin{tabular}{"+formato+"}\\toprule\n"
    texto += " & ".join(cabecera)+r" \\ \midrule"+"\n"
    texto += "\n".join(" & ".join(map(str, fila))+r" \\" for fila in filas)
    texto += "\n\\bottomrule\\end{tabular}\n"
    (SALIDA / "tablas" / nombre).write_text(texto, encoding="utf-8")


def guardar_figura(fig, nombre):
    fig.savefig(SALIDA/"figuras"/(nombre+".pdf"), metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(SALIDA/"figuras"/(nombre+".png"), dpi=170)
    plt.close(fig)


def generar():
    for carpeta in ("datos", "figuras", "tablas", "pdf"):
        (SALIDA/carpeta).mkdir(parents=True, exist_ok=True)
    x = np.linspace(0, X_MAX, 10001)
    referencia = aluminio(x)
    metricas = []
    for n in range(4):
        delta = diferencia(x, n)
        metricas.append([n, n+1, float(np.max(abs(delta))), float(np.sqrt(np.mean(delta**2))),
                         float(np.max(100*abs(np.expm1(delta)))),
                         float(np.max(abs(energia(x,n)-energia(x))))])
    escribir_csv("metricas.csv", ["orden", "terminos", "error_max_ln_gamma", "rmse_ln_gamma",
                                    "error_max_gamma_porcentaje", "error_max_gE_J_mol"], metricas)
    tabla_tex("metricas.tex", ["Orden", "Términos", "$E_{\\max}$", "RMSE", "$e_{\\gamma,\\max}$ (\\%)"],
              [[n,t,f"{e:.6f}",f"{r:.6f}",f"{g:.4f}"] for n,t,e,r,g,_ in metricas], "rrrrr")
    cabecera = ["x_Al"]
    columnas = [x]
    for n in range(4):
        p, h = aluminio(x,n), hierro(x,n)
        cabecera += [f"P{n}",f"gamma_Al_{n}",f"a_Al_{n}",f"H{n}",f"gE_{n}_J_mol"]
        columnas += [p,np.exp(p),x*np.exp(p),h,energia(x,n)]
    escribir_csv("malla.csv", cabecera, zip(*columnas))
    puntos = np.array([0.01,0.1,0.2,0.3,X_MAX])
    filas = list(zip(puntos, aluminio(puntos,1), aluminio(puntos,2), aluminio(puntos)))
    escribir_csv("composiciones.csv", ["x_Al","P1","P2","P3"], filas)
    tabla_tex("composiciones.tex", ["$x$","$P_1$","$P_2$","$F_3=P_3$"],
              [[f"{v:.6f}" for v in fila] for fila in filas], "rrrr")
    infinitas=[]
    for grado in (3,5,8,12):
        aproximacion=np.polynomial.polynomial.polyval(x,coeficientes_gamma(grado))
        infinitas.append([grado,float(np.max(100*abs(aproximacion/np.exp(referencia)-1)))])
    escribir_csv("taylor_gamma.csv", ["grado","error_max_porcentaje"], infinitas)
    tabla_tex("gamma.tex", ["Grado","Error relativo máximo (\\%)"],
              [[g,f"{e:.6f}"] for g,e in infinitas], "rr")
    resumen={"temperatura_K":T,"R_J_mol_K":R,"A":A,"B":B,"coeficientes_publicados":C.tolist(),
             "coeficientes_RK":C_RK.tolist(),"puntos":len(x),"x_max":X_MAX,
             "cruce_errores":C[2]/(-2*C[3]),
             "limites_uno_por_ciento":{str(n):limite_uno_por_ciento(n) for n in (1,2)},
             "gamma_Al_final":float(np.exp(aluminio(X_MAX))),
             "actividad_Al_final":float(X_MAX*np.exp(aluminio(X_MAX))),
             "ln_gamma_Fe_final":float(hierro(X_MAX)),
             "gE_final_kJ_mol":float(energia(X_MAX)/1000)}
    (SALIDA/"datos"/"resumen.json").write_text(json.dumps(resumen,indent=2,ensure_ascii=False)+"\n")
    comandos={"LimiteLineal":f"{limite_uno_por_ciento(1):.5f}",
              "LimiteCuadratico":f"{limite_uno_por_ciento(2):.5f}",
              "GammaFinal":f"{resumen['gamma_Al_final']:.6f}",
              "ActividadFinal":f"{resumen['actividad_Al_final']:.5f}",
              "HierroFinal":f"{resumen['ln_gamma_Fe_final']:.6f}",
              "EnergiaFinal":f"{resumen['gE_final_kJ_mol']:.5f}",
              "ErrorEnergiaLineal":f"{metricas[1][5]:.2f}",
              "ErrorEnergiaCuadratico":f"{metricas[2][5]:.2f}"}
    (SALIDA/"tablas"/"valores.tex").write_text("\n".join(
        "\\newcommand{\\"+k+"}{"+v+"}" for k,v in comandos.items())+"\n")
    from .desglose import generar_desglose
    generar_desglose(SALIDA)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,
                         "axes.spines.right":False,"axes.grid":True,"grid.alpha":0.2})
    colores=["#909090","#167d9a","#e89b2c","#753b83"]
    fig,axs=plt.subplots(1,2,figsize=(10,3.5),layout="constrained")
    for n in range(4):
        axs[0].plot(x,aluminio(x,n),label=f"Orden {n}",color=colores[n])
        axs[1].plot(x,np.exp(aluminio(x,n)),label=f"Orden {n}",color=colores[n])
    for ax,y in zip(axs,[r"$\ln\gamma_{Al}$",r"$\gamma_{Al}$"]):
        ax.set(xlabel=r"Fracción molar $x_{Al}$",ylabel=y); ax.legend(fontsize=8)
    guardar_figura(fig,"actividad")
    fig,axs=plt.subplots(1,2,figsize=(10,3.5),layout="constrained")
    xd=np.geomspace(1e-6,X_MAX,1000)
    for n in (0,1,2):
        axs[0].plot(x,100*abs(np.expm1(diferencia(x,n))),label=f"Orden {n}",color=colores[n])
        axs[1].loglog(xd,abs(diferencia(xd,n)),label=f"Orden {n}",color=colores[n])
    axs[0].axhline(1,color="black",ls=":",label="Tolerancia: 1 %")
    axs[0].set(ylabel="Error relativo en γAl (%)",ylim=(0,21))
    axs[1].set(ylabel="Error absoluto en ln γAl")
    for ax in axs:
        ax.set_xlabel(r"Fracción molar $x_{Al}$"); ax.legend(fontsize=8)
    guardar_figura(fig,"errores")
    fig,axs=plt.subplots(1,2,figsize=(10,3.5),layout="constrained")
    for n in (1,2,3):
        axs[0].plot(x,hierro(x,n),label=f"Orden {n}",color=colores[n])
        axs[1].plot(x,energia(x,n)/1000,label=f"Orden {n}",color=colores[n])
    for ax,y in zip(axs,[r"$\ln\gamma_{Fe}$",r"$g^E$ (kJ/mol)"]):
        ax.set(xlabel=r"Fracción molar $x_{Al}$",ylabel=y); ax.legend(fontsize=8)
    guardar_figura(fig,"termodinamica")


def compilar():
    if not shutil.which("pdflatex"):
        raise SystemExit("Falta pdflatex. Consulte README.md para instalar TeX Live o usar Docker.")
    entorno=os.environ.copy()
    entorno.update(SOURCE_DATE_EPOCH="1791331200", FORCE_SOURCE_DATE="1", TZ="UTC")
    for _ in range(3):
        resultado=subprocess.run(["pdflatex","-interaction=nonstopmode","-halt-on-error",
            "-output-directory=output/pdf","documentos/solucionario.tex"],
            cwd=RAIZ,env=entorno,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding="utf-8",errors="replace")
        if resultado.returncode:
            raise SystemExit(resultado.stdout)
    print("PDF generado: output/pdf/solucionario.pdf")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solo-datos",action="store_true",help="Generar cálculos y figuras sin compilar LaTeX")
    args=parser.parse_args()
    generar()
    if not args.solo_datos:
        compilar()
    archivos=sorted(p for p in SALIDA.rglob("*") if p.suffix in {".csv",".json",".png",".pdf",".tex"}
                    and p.name != "manifiesto.json" and "didactico" not in p.parts and "matlab" not in p.parts)
    manifiesto={str(p.relative_to(RAIZ)):hashlib.sha256(p.read_bytes()).hexdigest() for p in archivos}
    (SALIDA/"manifiesto.json").write_text(json.dumps(manifiesto,indent=2)+"\n")
    print("Resultados, figuras y huellas SHA-256 disponibles en output/.")
