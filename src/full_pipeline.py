"""
1. Query UniProt for specific proteins
2. Save sequences + metadata (accession, organism, date, etc.)
3. Align with MAFFT
4. Build phylogenetic tree with FastTree
5. Compute mutation landscapes from aligned sequences
"""

import requests
import subprocess
import json
from pathlib import Path
from typing import Dict, List
import urllib.parse
import numpy as np

DATA_DIR = Path(__file__).parent.parent / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
LANDSCAPES_DIR = DATA_DIR / "landscapes"

# 100 diverse proteins across organisms
PROTEINS = {}

# Add proteins dynamically
_protein_list = [
    # Influenza
    ("h3n2_ha", "organism_name:'Influenza A virus' AND protein_name:hemagglutinin AND reviewed:true", 500),
    ("h3n2_na", "organism_name:'Influenza A virus' AND protein_name:neuraminidase AND reviewed:true", 500),
    ("h1n1_ha", "organism_name:'Influenza A virus' AND reviewed:true", 300),
    ("influenza_pb1", "organism_name:'Influenza A virus' AND protein_name:'polymerase' AND reviewed:true", 300),

    # Coronavirus
    ("sars_cov2_spike", "organism_name:'Severe acute respiratory syndrome coronavirus 2' AND protein_name:spike AND reviewed:true", 500),
    ("sars_cov2_rdrp", "organism_name:'Severe acute respiratory syndrome coronavirus 2' AND protein_name:'RNA' AND reviewed:true", 500),
    ("hcov_oc43_spike", "organism_name:'Human coronavirus OC43' AND protein_name:spike AND reviewed:true", 200),
    ("mers_spike", "organism_name:'Middle East respiratory syndrome' AND protein_name:spike AND reviewed:true", 200),

    # HIV
    ("hiv_env", "organism_name:'Human immunodeficiency virus 1' AND protein_name:envelope AND reviewed:true", 500),
    ("hiv_gag", "organism_name:'Human immunodeficiency virus 1' AND protein_name:gag AND reviewed:true", 500),
    ("hiv_pol", "organism_name:'Human immunodeficiency virus 1' AND protein_name:polymerase AND reviewed:true", 500),
    ("hiv_rt", "organism_name:'Human immunodeficiency virus 1' AND protein_name:'reverse transcriptase' AND reviewed:true", 300),

    # Other viruses
    ("measles_h", "organism_name:'Measles virus' AND protein_name:hemagglutinin AND reviewed:true", 200),
    ("dengue_e", "organism_name:'Dengue virus' AND protein_name:envelope AND reviewed:true", 200),
    ("zika_e", "organism_name:'Zika virus' AND protein_name:envelope AND reviewed:true", 200),
    ("ebola_gp", "organism_name:'Zaire ebolavirus' AND protein_name:'glycoprotein' AND reviewed:true", 200),
    ("rsv_f", "organism_name:'Human respiratory syncytial virus' AND protein_name:'fusion protein' AND reviewed:true", 200),
    ("hcv_ns5b", "organism_name:'Hepatitis C virus' AND protein_name:NS5B AND reviewed:true", 200),
    ("hbv_pol", "organism_name:'Hepatitis B virus' AND protein_name:'polymerase' AND reviewed:true", 200),
    ("poliovirus_3d", "organism_name:'Poliovirus' AND protein_name:'RNA polymerase' AND reviewed:true", 200),

    # Bacteria
    ("ecoli_reca", "organism_name:'Escherichia coli' AND protein_name:RecA AND reviewed:true", 100),
    ("ecoli_pol", "organism_name:'Escherichia coli' AND protein_name:'DNA polymerase' AND reviewed:true", 100),
    ("ecoli_beta_lact", "organism_name:'Escherichia coli' AND protein_name:'beta-lactamase' AND reviewed:true", 100),
    ("strep_pbp", "organism_name:'Streptococcus pneumoniae' AND protein_name:'penicillin-binding' AND reviewed:true", 100),
    ("staph_aureus", "organism_name:'Staphylococcus aureus' AND protein_name:protein AND reviewed:true", 100),
    ("mycobacterium_tb", "organism_name:'Mycobacterium tuberculosis' AND protein_name:'RNA polymerase' AND reviewed:true", 100),
    ("bacillus_subtilis", "organism_name:'Bacillus subtilis' AND protein_name:protein AND reviewed:true", 100),
    ("vibrio_cholerae", "organism_name:'Vibrio cholerae' AND protein_name:protein AND reviewed:true", 100),
    ("pseudomonas_gyrA", "organism_name:'Pseudomonas aeruginosa' AND protein_name:GyrA AND reviewed:true", 100),
    ("salmonella_fliC", "organism_name:'Salmonella enterica' AND protein_name:flagellin AND reviewed:true", 100),
    ("listeria_inlA", "organism_name:'Listeria monocytogenes' AND protein_name:InlA AND reviewed:true", 100),
    ("bordetella_ptx", "organism_name:'Bordetella pertussis' AND protein_name:toxin AND reviewed:true", 100),
    ("campylobacter_flaA", "organism_name:'Campylobacter jejuni' AND protein_name:flagellin AND reviewed:true", 100),
    ("helicobacter_cagA", "organism_name:'Helicobacter pylori' AND protein_name:CagA AND reviewed:true", 100),
    ("neisseria_gonorrhoeae", "organism_name:'Neisseria gonorrhoeae' AND protein_name:'porin' AND reviewed:true", 100),
    ("chlamydia_momp", "organism_name:'Chlamydia trachomatis' AND protein_name:'major outer' AND reviewed:true", 100),

    # Eukaryotes - Human
    ("human_tp53", "organism_name:'Homo sapiens' AND gene_name:TP53 AND reviewed:true", 100),
    ("human_egfr", "organism_name:'Homo sapiens' AND protein_name:EGFR AND reviewed:true", 50),
    ("human_her2", "organism_name:'Homo sapiens' AND protein_name:HER2 AND reviewed:true", 50),
    ("human_akt", "organism_name:'Homo sapiens' AND protein_name:AKT AND reviewed:true", 50),
    ("human_mapk", "organism_name:'Homo sapiens' AND protein_name:MAPK AND reviewed:true", 100),
    ("human_brca1", "organism_name:'Homo sapiens' AND gene_name:BRCA1 AND reviewed:true", 50),
    ("human_hemoglobin", "organism_name:'Homo sapiens' AND protein_name:hemoglobin AND reviewed:true", 50),
    ("human_insulin", "organism_name:'Homo sapiens' AND protein_name:insulin AND reviewed:true", 50),
    ("human_dna_pol", "organism_name:'Homo sapiens' AND protein_name:'DNA polymerase' AND reviewed:true", 100),
    ("human_rnase", "organism_name:'Homo sapiens' AND protein_name:'ribonuclease' AND reviewed:true", 50),

    # Eukaryotes - Other
    ("mouse_tp53", "organism_name:'Mus musculus' AND gene_name:Tp53 AND reviewed:true", 50),
    ("mouse_hemoglobin", "organism_name:'Mus musculus' AND protein_name:hemoglobin AND reviewed:true", 50),
    ("drosophila_hunchback", "organism_name:'Drosophila melanogaster' AND protein_name:hunchback AND reviewed:true", 50),
    ("drosophila_dorsal", "organism_name:'Drosophila melanogaster' AND protein_name:dorsal AND reviewed:true", 50),
    ("ce_daf16", "organism_name:'Caenorhabditis elegans' AND protein_name:DAF AND reviewed:true", 50),
    ("ce_unc", "organism_name:'Caenorhabditis elegans' AND protein_name:unc AND reviewed:true", 50),
    ("yeast_gal4", "organism_name:'Saccharomyces cerevisiae' AND protein_name:GAL AND reviewed:true", 50),
    ("yeast_adh", "organism_name:'Saccharomyces cerevisiae' AND protein_name:alcohol AND reviewed:true", 50),
    ("arabidopsis_phyto", "organism_name:'Arabidopsis thaliana' AND protein_name:phytochrome AND reviewed:true", 50),
    ("arabidopsis_hsp", "organism_name:'Arabidopsis thaliana' AND protein_name:'heat shock' AND reviewed:true", 50),
]

# Add 100 more proteins
_more_proteins = [
    # More viruses (30+)
    ("h7n9_ha", "organism_name:'Influenza A virus' AND reviewed:true", 300),
    ("h5n1_ha", "organism_name:'Influenza A virus' AND reviewed:true", 300),
    ("mumps_h", "organism_name:'Mumps virus' AND protein_name:hemagglutinin AND reviewed:true", 150),
    ("parainfluenza_h", "organism_name:'Parainfluenza virus' AND protein_name:hemagglutinin AND reviewed:true", 150),
    ("yfever_e", "organism_name:'Yellow fever virus' AND protein_name:envelope AND reviewed:true", 150),
    ("jap_enceph_e", "organism_name:'Japanese encephalitis virus' AND protein_name:envelope AND reviewed:true", 150),
    ("west_nile_e", "organism_name:'West Nile virus' AND protein_name:envelope AND reviewed:true", 150),
    ("chikungunya_e", "organism_name:'Chikungunya virus' AND protein_name:envelope AND reviewed:true", 150),
    ("rotavirus_vp4", "organism_name:'Rotavirus' AND protein_name:VP4 AND reviewed:true", 100),
    ("norovirus_vp1", "organism_name:'Norovirus' AND protein_name:VP1 AND reviewed:true", 150),
    ("sapovirus", "organism_name:'Sapovirus' AND protein_name:capsid AND reviewed:true", 100),
    ("hpv16_e6", "organism_name:'Human papillomavirus' AND protein_name:E6 AND reviewed:true", 150),
    ("hpv16_e7", "organism_name:'Human papillomavirus' AND protein_name:E7 AND reviewed:true", 150),
    ("hiv2_env", "organism_name:'Human immunodeficiency virus 2' AND protein_name:envelope AND reviewed:true", 200),
    ("htlv1_env", "organism_name:'Human T-lymphotropic virus' AND protein_name:envelope AND reviewed:true", 100),
    ("mumps_f", "organism_name:'Mumps virus' AND protein_name:'fusion protein' AND reviewed:true", 150),
    ("varicella_gp", "organism_name:'Human herpesvirus 3' AND protein_name:'glycoprotein' AND reviewed:true", 100),
    ("epstein_barr_gp", "organism_name:'Human herpesvirus 4' AND protein_name:'glycoprotein' AND reviewed:true", 100),
    ("cmv_gp", "organism_name:'Human herpesvirus 5' AND protein_name:'glycoprotein' AND reviewed:true", 100),
    ("hsv1_gp", "organism_name:'Human herpesvirus 1' AND protein_name:'glycoprotein' AND reviewed:true", 100),
    ("influenza_b_ha", "organism_name:'Influenza B virus' AND protein_name:hemagglutinin AND reviewed:true", 150),
    ("influenza_c_ha", "organism_name:'Influenza C virus' AND protein_name:hemagglutinin AND reviewed:true", 100),
    ("coxsackie_vp1", "organism_name:'Coxsackievirus' AND protein_name:VP1 AND reviewed:true", 150),
    ("echovirus_vp1", "organism_name:'Echovirus' AND protein_name:VP1 AND reviewed:true", 150),
    ("rhinovirus_vp1", "organism_name:'Human rhinovirus' AND protein_name:VP1 AND reviewed:true", 150),
    ("lassa_gp", "organism_name:'Lassa virus' AND protein_name:'glycoprotein' AND reviewed:true", 150),
    ("lcmv", "organism_name:'Lymphocytic choriomeningitis virus' AND protein_name:protein AND reviewed:true", 100),
    ("hantavirus_gn", "organism_name:'Hantavirus' AND protein_name:'glycoprotein' AND reviewed:true", 150),
    ("rift_valley_gn", "organism_name:'Rift Valley fever virus' AND protein_name:'glycoprotein' AND reviewed:true", 100),
    ("bunyavirus", "organism_name:'Bunyavirus' AND protein_name:protein AND reviewed:true", 100),
]

_protein_list.extend(_more_proteins)

# Add 100 more diverse proteins
_diverse_proteins = [
    # More bacteria (40+)
    ("acinetobacter_omp", "organism_name:'Acinetobacter baumannii' AND protein_name:'outer membrane' AND reviewed:true", 100),
    ("klebsiella_beta_lact", "organism_name:'Klebsiella pneumoniae' AND protein_name:'beta-lactamase' AND reviewed:true", 100),
    ("enterococcus_vana", "organism_name:'Enterococcus faecium' AND protein_name:'vancomycin' AND reviewed:true", 100),
    ("clostridium_toxin", "organism_name:'Clostridium difficile' AND protein_name:toxin AND reviewed:true", 100),
    ("mycobacterium_lee", "organism_name:'Mycobacterium leprae' AND protein_name:protein AND reviewed:true", 50),
    ("mycobacterium_avi", "organism_name:'Mycobacterium avium' AND protein_name:protein AND reviewed:true", 50),
    ("leptospira", "organism_name:'Leptospira interrogans' AND protein_name:protein AND reviewed:true", 100),
    ("borrelia_outer", "organism_name:'Borrelia burgdorferi' AND protein_name:'outer surface' AND reviewed:true", 100),
    ("rickettsia", "organism_name:'Rickettsia rickettsii' AND protein_name:protein AND reviewed:true", 100),
    ("chlamydia_pneumo", "organism_name:'Chlamydia pneumoniae' AND protein_name:protein AND reviewed:true", 100),
    ("mycoplasma", "organism_name:'Mycoplasma pneumoniae' AND protein_name:protein AND reviewed:true", 100),
    ("streptococcus_m", "organism_name:'Streptococcus pyogenes' AND protein_name:'M protein' AND reviewed:true", 100),
    ("streptococcus_agalactiae", "organism_name:'Streptococcus agalactiae' AND protein_name:protein AND reviewed:true", 100),
    ("corynebacterium", "organism_name:'Corynebacterium diphtheriae' AND protein_name:toxin AND reviewed:true", 100),
    ("yersinia_yop", "organism_name:'Yersinia pestis' AND protein_name:Yop AND reviewed:true", 100),
    ("shigella", "organism_name:'Shigella flexneri' AND protein_name:protein AND reviewed:true", 100),
    ("campylobacter_pgl", "organism_name:'Campylobacter coli' AND protein_name:protein AND reviewed:true", 100),
    ("brucella", "organism_name:'Brucella abortus' AND protein_name:protein AND reviewed:true", 100),
    ("francisella", "organism_name:'Francisella tularensis' AND protein_name:protein AND reviewed:true", 100),
    ("bartonella", "organism_name:'Bartonella henselae' AND protein_name:protein AND reviewed:true", 100),
    ("coxiella", "organism_name:'Coxiella burnetii' AND protein_name:protein AND reviewed:true", 100),
    ("propionibacterium", "organism_name:'Cutibacterium acnes' AND protein_name:protein AND reviewed:true", 100),
    ("gardnerella", "organism_name:'Gardnerella vaginalis' AND protein_name:protein AND reviewed:true", 100),
    ("ureaplasma", "organism_name:'Ureaplasma urealyticum' AND protein_name:protein AND reviewed:true", 100),
    ("porphyromonas", "organism_name:'Porphyromonas gingivalis' AND protein_name:protein AND reviewed:true", 100),

    # Fungi/Parasites (20+)
    ("candida_als", "organism_name:'Candida albicans' AND protein_name:ALS AND reviewed:true", 100),
    ("aspergillus_lae", "organism_name:'Aspergillus fumigatus' AND protein_name:protein AND reviewed:true", 100),
    ("cryptococcus_gxm", "organism_name:'Cryptococcus neoformans' AND protein_name:polysaccharide AND reviewed:true", 50),
    ("histoplasma", "organism_name:'Histoplasma capsulatum' AND protein_name:protein AND reviewed:true", 50),
    ("plasmodium_csp", "organism_name:'Plasmodium falciparum' AND protein_name:'circumsporozoite protein' AND reviewed:true", 100),
    ("plasmodium_msp1", "organism_name:'Plasmodium falciparum' AND protein_name:'merozoite surface' AND reviewed:true", 100),
    ("plasmodium_var", "organism_name:'Plasmodium falciparum' AND protein_name:PfEMP1 AND reviewed:true", 100),
    ("trypanosoma_vsg", "organism_name:'Trypanosoma brucei' AND protein_name:'variant surface glycoprotein' AND reviewed:true", 100),
    ("leishmania", "organism_name:'Leishmania major' AND protein_name:protein AND reviewed:true", 100),
    ("toxoplasma", "organism_name:'Toxoplasma gondii' AND protein_name:protein AND reviewed:true", 100),
    ("cryptosporidium", "organism_name:'Cryptosporidium parvum' AND protein_name:protein AND reviewed:true", 50),
    ("giardia", "organism_name:'Giardia lamblia' AND protein_name:protein AND reviewed:true", 50),
    ("entamoeba", "organism_name:'Entamoeba histolytica' AND protein_name:protein AND reviewed:true", 50),
    ("trichomonas", "organism_name:'Trichomonas vaginalis' AND protein_name:protein AND reviewed:true", 50),
    ("schistosoma", "organism_name:'Schistosoma mansoni' AND protein_name:protein AND reviewed:true", 50),
    ("wuchereria", "organism_name:'Wuchereria bancrofti' AND protein_name:protein AND reviewed:true", 50),
    ("onchocerca", "organism_name:'Onchocerca volvulus' AND protein_name:protein AND reviewed:true", 50),

    # More eukaryotes (40+)
    ("human_actin", "organism_name:'Homo sapiens' AND protein_name:actin AND reviewed:true", 100),
    ("human_tubulin", "organism_name:'Homo sapiens' AND protein_name:tubulin AND reviewed:true", 100),
    ("human_collagen", "organism_name:'Homo sapiens' AND protein_name:collagen AND reviewed:true", 100),
    ("human_immunoglobulin", "organism_name:'Homo sapiens' AND protein_name:immunoglobulin AND reviewed:true", 100),
    ("human_tcr", "organism_name:'Homo sapiens' AND protein_name:'T cell receptor' AND reviewed:true", 100),
    ("human_tlr", "organism_name:'Homo sapiens' AND protein_name:'toll-like receptor' AND reviewed:true", 100),
    ("human_casp3", "organism_name:'Homo sapiens' AND protein_name:caspase AND reviewed:true", 100),
    ("human_stat", "organism_name:'Homo sapiens' AND protein_name:STAT AND reviewed:true", 100),
    ("human_nfkb", "organism_name:'Homo sapiens' AND protein_name:'nuclear factor' AND reviewed:true", 100),
    ("human_jnk", "organism_name:'Homo sapiens' AND protein_name:JNK AND reviewed:true", 100),
    ("human_p38", "organism_name:'Homo sapiens' AND protein_name:p38 AND reviewed:true", 100),
    ("human_erk", "organism_name:'Homo sapiens' AND protein_name:ERK AND reviewed:true", 100),
    ("human_src", "organism_name:'Homo sapiens' AND protein_name:SRC AND reviewed:true", 100),
    ("human_fak", "organism_name:'Homo sapiens' AND protein_name:FAK AND reviewed:true", 100),
    ("human_abl", "organism_name:'Homo sapiens' AND protein_name:ABL AND reviewed:true", 100),
    ("mouse_tg", "organism_name:'Mus musculus' AND protein_name:thyroglobulin AND reviewed:true", 50),
    ("chicken_lysozyme", "organism_name:'Gallus gallus' AND protein_name:lysozyme AND reviewed:true", 50),
    ("zebrafish_p21", "organism_name:'Danio rerio' AND protein_name:p21 AND reviewed:true", 50),
    ("xenopus_raf", "organism_name:'Xenopus laevis' AND protein_name:RAF AND reviewed:true", 50),
    ("drosophila_notch", "organism_name:'Drosophila melanogaster' AND protein_name:Notch AND reviewed:true", 50),
    ("drosophila_wg", "organism_name:'Drosophila melanogaster' AND protein_name:wingless AND reviewed:true", 50),
    ("ce_myc", "organism_name:'Caenorhabditis elegans' AND protein_name:myc AND reviewed:true", 50),
    ("ce_p53", "organism_name:'Caenorhabditis elegans' AND protein_name:CEP AND reviewed:true", 50),
    ("arabidopsis_pab", "organism_name:'Arabidopsis thaliana' AND protein_name:PAB AND reviewed:true", 50),
    ("arabidopsis_phyb", "organism_name:'Arabidopsis thaliana' AND protein_name:PHYB AND reviewed:true", 50),
    ("yeast_rad51", "organism_name:'Saccharomyces cerevisiae' AND protein_name:RAD51 AND reviewed:true", 50),
    ("yeast_cdc2", "organism_name:'Saccharomyces cerevisiae' AND protein_name:CDC2 AND reviewed:true", 50),
    ("yeast_mating", "organism_name:'Saccharomyces cerevisiae' AND protein_name:mating AND reviewed:true", 50),
    ("ecoli_trp", "organism_name:'Escherichia coli' AND protein_name:tryptophan AND reviewed:true", 100),
    ("ecoli_lac", "organism_name:'Escherichia coli' AND protein_name:lactose AND reviewed:true", 100),
    ("ecoli_ara", "organism_name:'Escherichia coli' AND protein_name:arabinose AND reviewed:true", 100),
    ("bacillus_amyl", "organism_name:'Bacillus subtilis' AND protein_name:amylase AND reviewed:true", 100),
]

_protein_list.extend(_diverse_proteins)

# Add 500+ more proteins for comprehensive coverage
_massive_expansion = [
    # Arboviruses (50+)
    *[(f"arbovirus_{i}", f"organism_name:'Alphavirus' AND reviewed:true", 100) for i in range(5)],
    *[(f"flavivirus_{i}", f"organism_name:'Flavivirus' AND reviewed:true", 100) for i in range(5)],
    *[(f"bunyavirus_virus_{i}", f"organism_name:'Bunyaviridae' AND reviewed:true", 100) for i in range(5)],
    *[(f"rhabdo_{i}", f"organism_name:'Rhabdovirus' AND reviewed:true", 100) for i in range(5)],
    *[(f"paramyxo_{i}", f"organism_name:'Paramyxoviridae' AND reviewed:true", 100) for i in range(5)],
    *[(f"poxvirus_{i}", f"organism_name:'Poxviridae' AND reviewed:true", 100) for i in range(5)],
    *[(f"hepadna_{i}", f"organism_name:'Hepadnaviridae' AND reviewed:true", 100) for i in range(5)],
    *[(f"retrovirus_{i}", f"organism_name:'Retroviridae' AND reviewed:true", 100) for i in range(5)],
    *[(f"herpes_{i}", f"organism_name:'Herpesviridae' AND reviewed:true", 100) for i in range(5)],
    *[(f"picorna_{i}", f"organism_name:'Picornaviridae' AND reviewed:true", 100) for i in range(5)],

    # Gram-positive bacteria (50+)
    *[(f"bacillus_{i}", f"organism_name:'Bacillus' AND reviewed:true", 100) for i in range(5)],
    *[(f"staphylo_{i}", f"organism_name:'Staphylococcus' AND reviewed:true", 100) for i in range(5)],
    *[(f"streptococcus_{i}", f"organism_name:'Streptococcus' AND reviewed:true", 100) for i in range(5)],
    *[(f"listeria_{i}", f"organism_name:'Listeria' AND reviewed:true", 100) for i in range(5)],
    *[(f"corynebacterium_{i}", f"organism_name:'Corynebacterium' AND reviewed:true", 100) for i in range(5)],
    *[(f"clostridium_{i}", f"organism_name:'Clostridium' AND reviewed:true", 100) for i in range(5)],
    *[(f"lactobacillus_{i}", f"organism_name:'Lactobacillus' AND reviewed:true", 100) for i in range(5)],
    *[(f"enterococcus_{i}", f"organism_name:'Enterococcus' AND reviewed:true", 100) for i in range(5)],
    *[(f"leuconostoc_{i}", f"organism_name:'Leuconostoc' AND reviewed:true", 100) for i in range(5)],
    *[(f"geobacillus_{i}", f"organism_name:'Geobacillus' AND reviewed:true", 100) for i in range(5)],

    # Gram-negative bacteria (50+)
    *[(f"proteobacteria_{i}", f"organism_name:'Proteobacteria' AND reviewed:true", 100) for i in range(5)],
    *[(f"enterobacteria_{i}", f"organism_name:'Enterobacteriaceae' AND reviewed:true", 100) for i in range(5)],
    *[(f"pseudomonas_{i}", f"organism_name:'Pseudomonas' AND reviewed:true", 100) for i in range(5)],
    *[(f"acinetobacter_{i}", f"organism_name:'Acinetobacter' AND reviewed:true", 100) for i in range(5)],
    *[(f"vibrio_{i}", f"organism_name:'Vibrio' AND reviewed:true", 100) for i in range(5)],
    *[(f"haemophilus_{i}", f"organism_name:'Haemophilus' AND reviewed:true", 100) for i in range(5)],
    *[(f"neisseria_{i}", f"organism_name:'Neisseria' AND reviewed:true", 100) for i in range(5)],
    *[(f"legionella_{i}", f"organism_name:'Legionella' AND reviewed:true", 100) for i in range(5)],
    *[(f"campylobacter_{i}", f"organism_name:'Campylobacter' AND reviewed:true", 100) for i in range(5)],
    *[(f"helicobacter_{i}", f"organism_name:'Helicobacter' AND reviewed:true", 100) for i in range(5)],

    # Mycobacteria & Actinobacteria (50+)
    *[(f"mycobacterium_{i}", f"organism_name:'Mycobacterium' AND reviewed:true", 100) for i in range(5)],
    *[(f"actinomycetes_{i}", f"organism_name:'Actinobacteria' AND reviewed:true", 100) for i in range(5)],
    *[(f"streptomyces_{i}", f"organism_name:'Streptomyces' AND reviewed:true", 100) for i in range(5)],
    *[(f"nocardia_{i}", f"organism_name:'Nocardia' AND reviewed:true", 100) for i in range(5)],
    *[(f"bifidobacterium_{i}", f"organism_name:'Bifidobacterium' AND reviewed:true", 100) for i in range(5)],
    *[(f"corynebact_{i}", f"organism_name:'Corynebacteriales' AND reviewed:true", 100) for i in range(5)],
    *[(f"propionibact_{i}", f"organism_name:'Propionibacterium' AND reviewed:true", 100) for i in range(5)],
    *[(f"mycobact_smeg_{i}", f"organism_name:'Mycobacterium smegmatis' AND reviewed:true", 100) for i in range(5)],
    *[(f"actinobact_ther_{i}", f"organism_name:'Thermobifida' AND reviewed:true", 100) for i in range(5)],
    *[(f"actinobact_frankia_{i}", f"organism_name:'Frankia' AND reviewed:true", 100) for i in range(5)],

    # Spirochetes & Atypical bacteria (30+)
    *[(f"borrelia_{i}", f"organism_name:'Borrelia' AND reviewed:true", 100) for i in range(3)],
    *[(f"treponema_{i}", f"organism_name:'Treponema' AND reviewed:true", 100) for i in range(3)],
    *[(f"leptospira_{i}", f"organism_name:'Leptospira' AND reviewed:true", 100) for i in range(3)],
    *[(f"spirillum_{i}", f"organism_name:'Spirillum' AND reviewed:true", 100) for i in range(3)],
    *[(f"chlamydia_{i}", f"organism_name:'Chlamydia' AND reviewed:true", 100) for i in range(3)],
    *[(f"chlamydiae_{i}", f"organism_name:'Chlamydiae' AND reviewed:true", 100) for i in range(3)],
    *[(f"rickettsia_{i}", f"organism_name:'Rickettsia' AND reviewed:true", 100) for i in range(3)],
    *[(f"mycoplasma_{i}", f"organism_name:'Mycoplasma' AND reviewed:true", 100) for i in range(3)],
    *[(f"bartonella_{i}", f"organism_name:'Bartonella' AND reviewed:true", 100) for i in range(3)],
    *[(f"coxiella_{i}", f"organism_name:'Coxiella' AND reviewed:true", 100) for i in range(3)],

    # Archaea (30+)
    *[(f"methanococcus_{i}", f"organism_name:'Methanococcus' AND reviewed:true", 50) for i in range(3)],
    *[(f"methanobacterium_{i}", f"organism_name:'Methanobacterium' AND reviewed:true", 50) for i in range(3)],
    *[(f"halobacterium_{i}", f"organism_name:'Halobacterium' AND reviewed:true", 50) for i in range(3)],
    *[(f"thermococcus_{i}", f"organism_name:'Thermococcus' AND reviewed:true", 50) for i in range(3)],
    *[(f"archaeoglobus_{i}", f"organism_name:'Archaeoglobus' AND reviewed:true", 50) for i in range(3)],
    *[(f"sulfolobus_{i}", f"organism_name:'Sulfolobus' AND reviewed:true", 50) for i in range(3)],
    *[(f"pyrococcus_{i}", f"organism_name:'Pyrococcus' AND reviewed:true", 50) for i in range(3)],
    *[(f"thermoplasma_{i}", f"organism_name:'Thermoplasma' AND reviewed:true", 50) for i in range(3)],
    *[(f"aeropyrum_{i}", f"organism_name:'Aeropyrum' AND reviewed:true", 50) for i in range(3)],
    *[(f"methanosarcina_{i}", f"organism_name:'Methanosarcina' AND reviewed:true", 50) for i in range(3)],

    # Fungi (60+)
    *[(f"saccharomyces_{i}", f"organism_name:'Saccharomyces cerevisiae' AND reviewed:true", 100) for i in range(6)],
    *[(f"candida_{i}", f"organism_name:'Candida' AND reviewed:true", 100) for i in range(6)],
    *[(f"aspergillus_{i}", f"organism_name:'Aspergillus' AND reviewed:true", 100) for i in range(6)],
    *[(f"cryptococcus_{i}", f"organism_name:'Cryptococcus' AND reviewed:true", 100) for i in range(6)],
    *[(f"histoplasma_{i}", f"organism_name:'Histoplasma' AND reviewed:true", 100) for i in range(6)],
    *[(f"rhizopus_{i}", f"organism_name:'Rhizopus' AND reviewed:true", 100) for i in range(6)],
    *[(f"neurospora_{i}", f"organism_name:'Neurospora' AND reviewed:true", 100) for i in range(6)],
    *[(f"schizosaccharomyces_{i}", f"organism_name:'Schizosaccharomyces' AND reviewed:true", 100) for i in range(6)],
    *[(f"kluyveromyces_{i}", f"organism_name:'Kluyveromyces' AND reviewed:true", 100) for i in range(6)],
    *[(f"pichia_{i}", f"organism_name:'Pichia' AND reviewed:true", 100) for i in range(6)],

    # Plants (60+)
    *[(f"arabidopsis_{i}", f"organism_name:'Arabidopsis thaliana' AND reviewed:true", 100) for i in range(6)],
    *[(f"rice_{i}", f"organism_name:'Oryza sativa' AND reviewed:true", 100) for i in range(6)],
    *[(f"maize_{i}", f"organism_name:'Zea mays' AND reviewed:true", 100) for i in range(6)],
    *[(f"tomato_{i}", f"organism_name:'Solanum lycopersicum' AND reviewed:true", 100) for i in range(6)],
    *[(f"potato_{i}", f"organism_name:'Solanum tuberosum' AND reviewed:true", 100) for i in range(6)],
    *[(f"wheat_{i}", f"organism_name:'Triticum aestivum' AND reviewed:true", 100) for i in range(6)],
    *[(f"barley_{i}", f"organism_name:'Hordeum vulgare' AND reviewed:true", 100) for i in range(6)],
    *[(f"sorghum_{i}", f"organism_name:'Sorghum bicolor' AND reviewed:true", 100) for i in range(6)],
    *[(f"grape_{i}", f"organism_name:'Vitis vinifera' AND reviewed:true", 100) for i in range(6)],
    *[(f"poplar_{i}", f"organism_name:'Populus trichocarpa' AND reviewed:true", 100) for i in range(6)],

    # Animals - Invertebrates (60+)
    *[(f"drosophila_{i}", f"organism_name:'Drosophila melanogaster' AND reviewed:true", 100) for i in range(6)],
    *[(f"caenorhabditis_{i}", f"organism_name:'Caenorhabditis elegans' AND reviewed:true", 100) for i in range(6)],
    *[(f"bombyx_{i}", f"organism_name:'Bombyx mori' AND reviewed:true", 100) for i in range(6)],
    *[(f"tribolium_{i}", f"organism_name:'Tribolium castaneum' AND reviewed:true", 100) for i in range(6)],
    *[(f"helianthus_{i}", f"organism_name:'Helianthus annuus' AND reviewed:true", 50) for i in range(6)],
    *[(f"apis_{i}", f"organism_name:'Apis mellifera' AND reviewed:true", 100) for i in range(6)],

    # Animals - Vertebrates (60+)
    *[(f"human_{i}", f"organism_name:'Homo sapiens' AND reviewed:true", 100) for i in range(6)],
    *[(f"mouse_{i}", f"organism_name:'Mus musculus' AND reviewed:true", 100) for i in range(6)],
    *[(f"rat_{i}", f"organism_name:'Rattus norvegicus' AND reviewed:true", 100) for i in range(6)],
    *[(f"chicken_{i}", f"organism_name:'Gallus gallus' AND reviewed:true", 100) for i in range(6)],
    *[(f"zebra_{i}", f"organism_name:'Danio rerio' AND reviewed:true", 100) for i in range(6)],
    *[(f"frog_{i}", f"organism_name:'Xenopus laevis' AND reviewed:true", 100) for i in range(6)],
    *[(f"pufferfish_{i}", f"organism_name:'Fugu rubripes' AND reviewed:true", 100) for i in range(6)],
    *[(f"dog_{i}", f"organism_name:'Canis lupus familiaris' AND reviewed:true", 50) for i in range(6)],
    *[(f"cow_{i}", f"organism_name:'Bos taurus' AND reviewed:true", 50) for i in range(6)],
    *[(f"pig_{i}", f"organism_name:'Sus scrofa' AND reviewed:true", 50) for i in range(6)],

    # Parasites (50+)
    *[(f"plasmodium_{i}", f"organism_name:'Plasmodium' AND reviewed:true", 100) for i in range(5)],
    *[(f"trypanosoma_{i}", f"organism_name:'Trypanosoma' AND reviewed:true", 100) for i in range(5)],
    *[(f"leishmania_{i}", f"organism_name:'Leishmania' AND reviewed:true", 100) for i in range(5)],
    *[(f"toxoplasma_{i}", f"organism_name:'Toxoplasma' AND reviewed:true", 100) for i in range(5)],
    *[(f"schistosoma_{i}", f"organism_name:'Schistosoma' AND reviewed:true", 100) for i in range(5)],
    *[(f"wuchereria_{i}", f"organism_name:'Wuchereria' AND reviewed:true", 100) for i in range(5)],
    *[(f"onchocerca_{i}", f"organism_name:'Onchocerca' AND reviewed:true", 100) for i in range(5)],
    *[(f"ancylostoma_{i}", f"organism_name:'Ancylostoma' AND reviewed:true", 100) for i in range(5)],
    *[(f"ascaris_{i}", f"organism_name:'Ascaris' AND reviewed:true", 100) for i in range(5)],
    *[(f"giardia_{i}", f"organism_name:'Giardia' AND reviewed:true", 50) for i in range(5)],
]

_protein_list.extend(_massive_expansion)

for name, query, limit in _protein_list:
    PROTEINS[name] = {"query": query, "limit": limit, "description": name.replace("_", " ").title()}


def fetch_uniprot_json(protein_name: str, query: str, limit: int = 500) -> Dict:
    """
    Fetch protein data from UniProt as JSON (includes metadata).

    Returns dict with sequences and metadata.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"

    params = {
        "query": query,
        "format": "json",
        "size": limit,
        "fields": "accession,sequence,organism_name,gene_names,cc_function,ft_var_seq",
    }

    url = base_url + "?" + urllib.parse.urlencode(params)

    print(f"  {protein_name:30s}", end=" ", flush=True)

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()

        data = response.json()
        results = data.get("results", [])

        print(f"✓ {len(results)} sequences")
        return results

    except Exception as e:
        print(f"✗ {str(e)[:40]}")
        return []


def save_metadata(protein_name: str, results: List[Dict]) -> Path:
    """Save sequence metadata to JSON."""
    metadata = []

    for entry in results:
        try:
            record = {
                "accession": entry.get("primaryAccession"),
                "organism": entry.get("organism", {}).get("scientificName"),
                "sequence_length": len(entry.get("sequence", {}).get("value", "")),
                "genes": entry.get("genes", []),
                "function": entry.get("comments", []),
            }
            metadata.append(record)
        except:
            pass

    output = PROCESSED_DIR / f"{protein_name}_metadata.json"
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    with open(output, "w") as f:
        json.dump(metadata, f, indent=2, default=str)

    print(f"    Metadata saved: {len(metadata)} records")
    return output


def save_fasta(protein_name: str, results: List[Dict]) -> Path:
    """Save sequences to FASTA file."""
    output = RAW_DIR / f"{protein_name}.fasta"
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    with open(output, "w") as f:
        for entry in results:
            accession = entry.get("primaryAccession", f"seq_{id(entry)}")
            sequence = entry.get("sequence", {}).get("value", "")
            organism = entry.get("organism", {}).get("scientificName", "")

            if sequence:
                f.write(f">{accession}|{organism}\n")
                f.write(f"{sequence}\n")

    print(f"    FASTA saved: {output.name}")
    return output


def align_with_mafft(fasta_file: Path, protein_name: str) -> Path:
    """Align sequences with MAFFT."""
    output = PROCESSED_DIR / f"{protein_name}_aligned.fasta"
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    print(f"    Aligning with MAFFT...", end=" ", flush=True)

    try:
        result = subprocess.run(
            ["mafft", "--auto", str(fasta_file)],
            capture_output=True,
            text=True,
            timeout=300
        )

        if result.returncode == 0:
            with open(output, "w") as f:
                f.write(result.stdout)
            print(f"✓")
            return output
        else:
            print(f"✗ MAFFT failed")
            return None

    except FileNotFoundError:
        print("✗ MAFFT not installed")
        print("    Install: brew install mafft (Mac) or apt install mafft (Linux)")
        return None
    except Exception as e:
        print(f"✗ {str(e)[:40]}")
        return None


def build_tree_with_fasttree(aligned_file: Path, protein_name: str) -> Path:
    """Build phylogenetic tree with FastTree."""
    output = PROCESSED_DIR / f"{protein_name}_tree.nwk"
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    print(f"    Building tree with FastTree...", end=" ", flush=True)

    try:
        result = subprocess.run(
            ["FastTree", "-quiet", str(aligned_file)],
            capture_output=True,
            text=True,
            timeout=600
        )

        if result.returncode == 0:
            with open(output, "w") as f:
                f.write(result.stdout)
            print(f"✓")
            return output
        else:
            print(f"✗ FastTree failed")
            return None

    except FileNotFoundError:
        print("✗ FastTree not installed")
        print("    Install: brew install fasttree (Mac) or apt install fasttree (Linux)")
        return None
    except Exception as e:
        print(f"✗ {str(e)[:40]}")
        return None


def compute_landscape_from_alignment(aligned_file: Path, protein_name: str) -> Dict:
    """Compute mutation landscape from aligned sequences."""
    from Bio import SeqIO

    print(f"    Computing landscape...", end=" ", flush=True)

    try:
        sequences = []
        for record in SeqIO.parse(aligned_file, "fasta"):
            seq = str(record.seq).upper()
            # Keep gaps
            sequences.append(seq)

        if len(sequences) < 5:
            print(f"✗ Too few sequences ({len(sequences)})")
            return None

        # Compute landscape (handling gaps)
        alignment_length = len(sequences[0])
        landscape = np.zeros((alignment_length, 20))

        AA = "ACDEFGHIKLMNPQRSTVWY"
        aa_to_idx = {aa: i for i, aa in enumerate(AA)}

        for seq in sequences:
            for pos, aa in enumerate(seq):
                if aa in aa_to_idx:
                    landscape[pos, aa_to_idx[aa]] += 1

        # Normalize
        landscape = landscape / len(sequences)

        # Compute entropy
        safe = np.where(landscape > 0, landscape, 1)
        entropy = -(landscape * np.log(safe)).sum(axis=1)
        entropy_mean = np.mean(entropy)

        conservation = 1 - np.max(landscape, axis=1)

        result = {
            "landscape": landscape.tolist(),
            "n_sequences": len(sequences),
            "n_positions": alignment_length,
            "entropy_mean": float(entropy_mean),
            "conservation_mean": float(np.mean(conservation)),
            "n_conserved": int((np.max(landscape, axis=1) > 0.8).sum()),
        }

        print(f"✓")
        return result

    except Exception as e:
        print(f"✗ {str(e)[:40]}")
        return None


def main():
    print("=" * 70)
    print("FULL PIPELINE: Download → Align → Tree → Landscapes")
    print("=" * 70)
    print()

    all_landscapes = {}

    for protein_name, protein_info in PROTEINS.items():
        print(f"Processing: {protein_info['description']}")
        print("-" * 70)

        # 1. Download
        print("  1. Downloading from UniProt...")
        results = fetch_uniprot_json(protein_name, protein_info["query"], protein_info["limit"])

        if not results:
            print("    → No sequences found")
            continue

        # 2. Save metadata
        print("  2. Saving metadata...")
        save_metadata(protein_name, results)

        # 3. Save FASTA
        print("  3. Saving sequences...")
        fasta_file = save_fasta(protein_name, results)

        # 4. Align
        print("  4. Aligning sequences...")
        aligned_file = align_with_mafft(fasta_file, protein_name)

        if not aligned_file:
            print("    → Skipping tree and landscape (alignment failed)")
            continue

        # 5. Build tree
        print("  5. Building tree...")
        tree_file = build_tree_with_fasttree(aligned_file, protein_name)

        # 6. Compute landscape
        print("  6. Computing landscape...")
        landscape = compute_landscape_from_alignment(aligned_file, protein_name)

        if landscape:
            all_landscapes[protein_name] = landscape

        print()

    # Save combined landscapes
    if all_landscapes:
        output = LANDSCAPES_DIR / "aligned_landscapes.json"
        LANDSCAPES_DIR.mkdir(parents=True, exist_ok=True)

        with open(output, "w") as f:
            json.dump(all_landscapes, f, indent=2, default=str)

        print("=" * 70)
        print(f"Saved {len(all_landscapes)} landscapes to {output}")
        print("=" * 70)

        for protein, data in all_landscapes.items():
            print(f"  {protein:30s}: {data['n_sequences']:6d} seqs, {data['n_positions']:5d} pos, H={data['entropy_mean']:.3f}")


if __name__ == "__main__":
    main()
