# 04. Phylogenetic signal detection: LQ

Analisi riproducibile del segnale filogenetico del longevity quotient non
trasformato (`LQ_mammal`). Calcola K di Blomberg con test di randomizzazione e
λ di Pagel con likelihood-ratio test rispetto a λ = 0.

La cartella è autosufficiente: contiene copie locali del fenotipo LQ prodotto
dall'analisi 02 e dell'albero Kuderna S4 usato nell'analisi 03.

## Esecuzione

Sono richiesti R e i pacchetti `ape` e `phytools`. Lo script risolve gli input
rispetto alla propria posizione e può quindi essere lanciato da qualunque
directory:

```bash
Rscript /percorso/04.phylognetic.signal.detection/run_phylogenetic_signal.R
```

Il file `results.md` viene rigenerato a ogni esecuzione. Il seed del test di
randomizzazione, il numero di permutazioni, le versioni dei pacchetti e le
impronte degli input sono registrati nel report.
