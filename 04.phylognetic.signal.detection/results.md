# 04. Phylogenetic signal detection: LQ

Report generato il 2026-08-17.

## Risultati

| Statistica | Stima | Test dell'assenza di segnale | P-value | Significativo (α = 0.05) |
|---|---:|---|---:|:---:|
| K di Blomberg | 0.145606 | Randomizzazione delle etichette (9999 permutazioni) | 3.0003e-04 | sì |
| λ di Pagel | 0.785893 | Likelihood-ratio test, H0: λ = 0 (LR = 34.6521) | 3.9422e-09 | sì |

## Interpretazione

- K è 0.145606: inferiore a 1, quindi la somiglianza tra parenti è più debole di quella attesa sotto moto browniano. Il test di randomizzazione rifiuta l'assenza di segnale filogenetico.
- λ è 0.785893; 0 indica assenza di covarianza filogenetica e valori vicini a 1 una struttura simile all'attesa browniana. Il likelihood-ratio test rifiuta λ = 0.
- **Conclusione:** Entrambi i test rilevano un segnale filogenetico significativo nel LQ.

## Dati e controlli

- Fenotipo: `LQ_mammal` (LQ non trasformato).
- Specie nel dataset: 162.
- Specie nell'albero completo: 236.
- Specie analizzate dopo l'intersezione: 126.
- Specie del dataset escluse perché assenti dall'albero: 36.
- Punte dell'albero escluse perché prive di LQ: 110.
- Albero analizzato: radicato = sì; biforcante = sì; ultrametrico = sì.

Le specie senza corrispondenza esatta nell'albero sono:

`Aotus_lemurinus`, `Ateles_fusciceps`, `Brachyteles_arachnoides`, `Callithrix_penicillata`, `Cebus_capucinus`, `Cercocebus_agilis`, `Cercocebus_atys`, `Cercopithecus_campbelli`, `Cercopithecus_erythrogaster`, `Cercopithecus_wolfi`, `Chiropotes_satanas`, `Chlorocebus_aethiops`, `Galagoides_demidovii`, `Galagoides_zanzibaricus`, `Hylobates_moloch`, `Hylobates_pileatus`, `Lagothrix_lagotricha`, `Leontopithecus_chrysopygus`, `Lepilemur_mustelinus`, `Lophocebus_albigena`, `Macaca_ochreata`, `Macaca_sinica`, `Macaca_sylvanus`, `Miopithecus_talapoin`, `Mirza_coquereli`, `Nomascus_leucogenys`, `Phaner_furcifer`, `Pithecia_monachus`, `Pithecia_pithecia`, `Plecturocebus_donacophilus`, `Presbytis_melalophos`, `Saguinus_leucopus`, `Saimiri_boliviensis`, `Symphalangus_syndactylus`, `Tarsius_bancanus`, `Tarsius_tarsier`.

## Metodo

I nomi delle specie sono stati abbinati esattamente tra `accepted_tree_tip` e le punte dell'albero. L'albero è stato potato alle specie con LQ disponibile. K di Blomberg è stato testato permutando i valori del tratto sulle punte; λ di Pagel è stato stimato per massima verosimiglianza e confrontato con λ = 0 mediante likelihood-ratio test. Entrambe le statistiche sono state calcolate con `phytools::phylosig`.

## Riproducibilità

Dalla cartella dell'analisi:

```bash
Rscript run_phylogenetic_signal.R
```

- Seed del test di randomizzazione: `20260817`.
- R: `R version 4.3.1 Patched (2023-08-12 r84939)`.
- ape: `5.8.1`; phytools: `2.5.2`.
- MD5 `input/longevity_quotient.tsv`: `69f3a39097554167716c5c5002220ce0`.
- MD5 `input/science.abn7829_data_s4.nex.tree`: `e6793ca5602c702b81ca255d329ee01b`.

Gli input sono copie locali incluse nella cartella, quindi lo script non dipende da file esterni al modulo di analisi.
