import os
import csv

# Définition des chemins relatifs
base_dir = "results/PRODUCTION_HYBRID"
out_file = "results/PRODUCTION_HYBRID/matam_salmon_master_export.tsv"

print(f"Création du fichier unique : {out_file}...")

# Ouverture du fichier de sortie en mode écriture ("w"). 
# L'argument newline="" est une bonne pratique avec la librairie csv en Python 3, 
# il évite l'apparition de sauts de lignes parasites entre les entrées sous certains OS.
with open(out_file, "w", newline="") as f_out:
    writer = csv.writer(f_out, delimiter="\t")
    
    # Création des en-têtes (Header)
    # Ce format "long" est parfait pour R (tidy data). Chaque ligne représentera 
    # l'observation d'un taxon spécifique, dans un échantillon spécifique, avec sa séquence.
    writer.writerow(["Sample", "Scaffold_ID", "NumReads", "Sequence"])

    # Parcours itératif et trié des dossiers d'échantillons.
    # Le tri (sorted) garantit que les échantillons apparaîtront dans le même ordre à chaque exécution.
    for sample in sorted(os.listdir(base_dir)):
        sample_dir = os.path.join(base_dir, sample)
        
        # Sécurité : on ignore silencieusement les fichiers isolés (comme des logs) 
        # qui pourraient traîner à la racine du dossier PRODUCTION_HYBRID.
        if not os.path.isdir(sample_dir): 
            continue

        # Reconstruction dynamique des chemins vers les livrables de l'échantillon courant
        fasta_path = os.path.join(sample_dir, f"{sample}_scaffolds.fasta")
        salmon_path = os.path.join(sample_dir, f"{sample}_abundance_salmon.tsv")

        # Vérification de l'intégrité de l'analyse : on ne procède à la fusion que 
        # si les deux étapes précédentes (MATAM et Salmon) ont généré leurs fichiers.
        if os.path.exists(fasta_path) and os.path.exists(salmon_path):
            
            # =========================================================================
            # 1. Parsing du fichier FASTA (Mise en cache des séquences)
            # =========================================================================
            # On utilise un dictionnaire car il offre une complexité de recherche O(1).
            # C'est parfait pour associer rapidement un ID de scaffold à sa séquence plus tard.
            seqs = {}
            current_id = ""
            with open(fasta_path, "r") as fa:
                for line in fa:
                    line = line.strip() # Nettoyage des retours chariots (\n) et espaces
                    if line.startswith(">"):
                        # Extraction de l'identifiant pur (ex: ">scaffold_1 length=1500" devient "scaffold_1")
                        current_id = line[1:].split()[0]
                        seqs[current_id] = []
                    else:
                        # Stockage progressif des nucléotides. Utiliser .append() sur une liste 
                        # est algorithmiquement plus efficace en Python que de concaténer des chaînes (+=) 
                        # ligne par ligne, surtout pour de longs scaffolds 16S.
                        seqs[current_id].append(line)
                        
            # Reconstitution des séquences complètes une fois le fichier lu
            for k in seqs: 
                seqs[k] = "".join(seqs[k])

            # =========================================================================
            # 2. Parsing du fichier Salmon (Abondances) et Écriture fusionnée
            # =========================================================================
            with open(salmon_path, "r") as sf:
                # csv.DictReader est le choix optimal ici : il lit automatiquement la première 
                # ligne comme en-tête et permet d'extraire les données par le nom de la colonne 
                # ("Name", "NumReads") plutôt que par un index abstrait (row[0], row[4]).
                reader = csv.DictReader(sf, delimiter="\t")
                for row in reader:
                    scaf_id = row["Name"]
                    
                    # NumReads est casté en 'float'. Pourquoi ? Salmon utilise un algorithme 
                    # Expectation-Maximization probabiliste. Il peut estimer qu'un read appartient 
                    # à 70% au taxon A et 30% au taxon B, générant des abondances non-entières (ex: 15.4 reads).
                    reads = float(row["NumReads"])
                    
                    # Récupération de la séquence associée via le dictionnaire. 
                    # .get() renvoie une chaîne vide "" si l'ID n'est pas trouvé (sécurité).
                    sequence = seqs.get(scaf_id, "")
                    
                    # Filtre biologique et optimisation : on ignore silencieusement les scaffolds 
                    # ayant une abondance de 0 (ou absents du FASTA). Cela réduit drastiquement 
                    # le "bruit" et la taille du fichier TSV final.
                    if sequence and reads > 0:
                        writer.writerow([sample, scaf_id, reads, sequence])
                        
            print(f"  -> {sample} ajouté.")

print("Terminé avec succès !")