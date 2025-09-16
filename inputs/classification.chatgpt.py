from collections import defaultdict
import pandas as pd

# Lista dei generi fornita
genera = [
    "Allenopithecus", "Allochrocebus", "Alouatta", "Aotus", "Arctocebus", "Ateles", "Avahi", 
    "Brachyteles", "Cacajao", "Callibella", "Callicebus", "Callimico", "Callithrix", "Carlito", 
    "Cebuella", "Cebus", "Cephalopachus", "Cercocebus", "Cercopithecus", "Cheirogaleus", 
    "Cheracebus", "Chiropotes", "Chlorocebus", "Colobus", "Daubentonia", "Erythrocebus", 
    "Eulemur", "Galago", "Galagoides", "Gorilla", "Hapalemur", "Hoolock", "Hylobates", "Indri", 
    "Lagothrix", "Lemur", "Leontocebus", "Leontopithecus", "Lepilemur", "Lophocebus", "Loris", 
    "Macaca", "Mandrillus", "Mico", "Microcebus", "Miopithecus", "Mirza", "Nasalis", "Nomascus", 
    "Nycticebus", "Otolemur", "Pan", "Papio", "Paragalago", "Perodicticus", "Piliocolobus", 
    "Pithecia", "Plecturocebus", "Pongo", "Presbytis", "Prolemur", "Propithecus", "Ptilocercus", 
    "Pygathrix", "Rhinopithecus", "Saguinus", "Saimiri", "Sapajus", "Semnopithecus", 
    "Symphalangus", "Tarsius", "Theropithecus", "Trachypithecus", "Unknown", "Varecia"
]

# Famiglie tassonomiche note per ciascun genere (manual mapping)
family_map = {
    # Strepsirrhini
    "Cheirogaleus": "Cheirogaleidae", "Microcebus": "Cheirogaleidae", "Mirza": "Cheirogaleidae", "Allocebus": "Cheirogaleidae",
    "Lepilemur": "Lepilemuridae",
    "Avahi": "Indriidae", "Propithecus": "Indriidae", "Indri": "Indriidae",
    "Daubentonia": "Daubentoniidae",
    "Lemur": "Lemuridae", "Varecia": "Lemuridae", "Eulemur": "Lemuridae", "Hapalemur": "Lemuridae", "Prolemur": "Lemuridae",
    "Galago": "Galagidae", "Galagoides": "Galagidae", "Otolemur": "Galagidae", "Paragalago": "Galagidae",
    "Loris": "Lorisidae", "Nycticebus": "Lorisidae", "Perodicticus": "Lorisidae", "Arctocebus": "Lorisidae",
    "Tarsius": "Tarsiidae", "Cephalopachus": "Tarsiidae",  # Tarsiers
    "Ptilocercus": "Ptilocercidae",  # Not a primate but part of tree shrews

    # Platyrrhini
    "Callithrix": "Callitrichidae", "Mico": "Callitrichidae", "Callimico": "Callitrichidae", "Saguinus": "Callitrichidae", "Leontopithecus": "Callitrichidae", "Leontocebus": "Callitrichidae", "Cebuella": "Callitrichidae", "Callibella": "Callitrichidae",
    "Aotus": "Aotidae",
    "Saimiri": "Cebidae", "Cebus": "Cebidae", "Sapajus": "Cebidae",
    "Callicebus": "Pitheciidae", "Plecturocebus": "Pitheciidae", "Cheracebus": "Pitheciidae",
    "Pithecia": "Pitheciidae", "Chiropotes": "Pitheciidae", "Cacajao": "Pitheciidae",
    "Alouatta": "Atelidae", "Ateles": "Atelidae", "Brachyteles": "Atelidae", "Lagothrix": "Atelidae",

    # Catarrhini
    "Macaca": "Cercopithecidae", "Papio": "Cercopithecidae", "Mandrillus": "Cercopithecidae", "Theropithecus": "Cercopithecidae",
    "Cercopithecus": "Cercopithecidae", "Chlorocebus": "Cercopithecidae", "Erythrocebus": "Cercopithecidae", "Miopithecus": "Cercopithecidae",
    "Allenopithecus": "Cercopithecidae", "Lophocebus": "Cercopithecidae", "Cercocebus": "Cercopithecidae",
    "Colobus": "Cercopithecidae", "Piliocolobus": "Cercopithecidae", "Nasalis": "Cercopithecidae", "Simias": "Cercopithecidae",
    "Presbytis": "Cercopithecidae", "Pygathrix": "Cercopithecidae", "Rhinopithecus": "Cercopithecidae", "Semnopithecus": "Cercopithecidae", "Trachypithecus": "Cercopithecidae",

    "Hylobates": "Hylobatidae", "Nomascus": "Hylobatidae", "Symphalangus": "Hylobatidae", "Hoolock": "Hylobatidae",
    "Pan": "Hominidae", "Gorilla": "Hominidae", "Pongo": "Hominidae", "Homo": "Hominidae",
    "Carlito": "Tarsiidae",  # sometimes listed separately
    "Unknown": "Unknown"
}

# Raccogli in un dizionario famiglia -> lista generi
family_dict = defaultdict(list)
for genus in genera:
    family = family_map.get(genus, "Unknown")
    family_dict[family].append(genus)

# Prepara il tabulato
df = pd.DataFrame({
    "family": list(family_dict.keys()),
    "genera": [",".join(sorted(glist)) for glist in family_dict.values()]
})

import ace_tools as tools; tools.display_dataframe_to_user(name="Genera per Famiglia", dataframe=df)


