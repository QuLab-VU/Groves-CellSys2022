arcCols <- c("lipogenesis"="#ffff00",
	     "biomass & energy"="#019e59", 
	     "cell division"="#007eb1", 
	     "immune interaction"="#d02690", 
	     "invasion & signaling"="#111111")

## List of commonly enriched MSigDB pathways, representative of each cancer hallmark
toShow <-
    list(
        "Resisting cell death"=c(
            "ANTI APOPTOSIS",
            "REACTOME EXTENSION OF TELOMERES"
        ),
        "Evading growth suppressors"=c(
            "REACTOME DNA REPLICATION",
            "M PHASE OF MITOTIC CELL CYCLE",
            "KEGG CELL CYCLE"
        ),
        "Interaction with immune system"=c(
            "BIOCARTA INFLAM PATHWAY",
            "BIOCARTA TCYTOTOXIC PATHWAY",
            "REACTOME INTERFERON GAMMA SIGNALING",
            "REACTOME PD1 SIGNALING",
            "BIOCARTA CTLA4 PATHWAY",
            "KEGG ALLOGRAFT REJECTION",
            "BIOCARTA LECTIN PATHWAY"),
        "Energetics"=c(
            "REACTOME GLYCOLYSIS",
            "REACTOME RESPIRATORY ELECTRON TRANSPORT",
            "KEGG RIBOSOME",
            "MITOCHONDRIAL RIBOSOME",
            "KEGG PROTEASOME"
	),
        "DNA repair"=c(
            "REACTOME DOUBLE STRAND BREAK REPAIR",
            "KEGG NUCLEOTIDE EXCISION REPAIR"
        ),
        "Angiogenesis"=c(
            "PID VEGF VEGFR PATHWAY",
            "POSITIVE REGULATION OF ANGIOGENESIS"),
        "Signaling"=c(
            "PID HEDGEHOG 2PATHWAY",
            "INSULIN LIKE GROWTH FACTOR RECEPTOR BINDING",
            "REACTOME FGFR LIGAND BINDING AND ACTIVATION",
            "PID WNT SIGNALING PATHWAY",
            "KEGG TGF BETA SIGNALING PATHWAY",
            "PID RAS PATHWAY",
            "BIOCARTA MTOR PATHWAY"
        ),
        "Peroxisome"=c(
            "REACTOME PEROXISOMAL LIPID METABOLISM",
            "KEGG GLYCOSYLPHOSPHATIDYLINOSITOL GPI ANCHOR BIOSYNTHESIS",
            "KEGG ASCORBATE AND ALDARATE METABOLISM"
        ),
        "Invasion & metastasis"=c(
            "REACTOME DEGRADATION OF THE EXTRACELLULAR MATRIX",
            "CELL MIGRATION"
        ),
        "Micro-environment"=c(
            "TRANSMISSION OF NERVE IMPULSE",
            "REGULATION OF NEUROGENESIS",
            "COLLAGEN",
            "REACTOME SMOOTH MUSCLE CONTRACTION",
	    "RESPONSE TO HYPOXIA"
        )
    )

fmtMSigDBnames <- function(MSigDBnames, maxWidth=35) {
    MSigDBnames %>% 
        str_replace("^REACTOME ", "") %>%
        str_replace("^KEGG ", "") %>%
        str_replace("^PID ", "") %>%
        str_replace("^BIOCARTA ", "") %>%
        str_replace("THE ", " ") %>%
        str_to_title() %>%
        str_trunc(width=maxWidth)
}

lipidShow <-
    list(
        "Fatty acid synthesis"=c(
	    "SREBF1", #TF controlling FASN expression
            "ACLY", #cytosolic ATP cytrate lyase
            "ACACA", #Acetyl-CoA carboxylase
            "FASN", #fatty acid synthase
            "SCD" #stearoyl-CoA desaturase
        ),
	"anaerobic glycosis"=c(
	    "LDHA", #lactate dehydrogenase
	    "LDHB" #lactate dehydrogenase
	),
	"NADPH generation"=c(
	    "ME1", #malic enzyme
	    "ME2" #malic enzyme
	),
        "beta-oxidation"=c(
            "ACOX1",
            "CPT1B"
        ),
        "key peroxisome proteins"=c(
            "CAT",
            "PEX19"
        ),
        "other"=c(
            "GPX4"
        )
    )

##################################################

## Shorter list of MSigDB pathways for each cancer hallmark

toShowFig <-
    list(
        "Resisting cell death"=c(
            "ANTI APOPTOSIS"
        ),
        "Sustaining proliferative signaling"=c(
            "REACTOME PI3K CASCADE",
            "PID RAS PATHWAY"
        ),
        "Evading growth suppressors"=c(
            "REACTOME DNA REPLICATION",
            "MITOSIS"
        ),
        "Activating invasion & metastasis"=c(
            "REACTOME DEGRADATION OF THE EXTRACELLULAR MATRIX",
            "CELL MIGRATION"
        ),
        "Enabling replicative immortability"=c(
            "REACTOME EXTENSION OF TELOMERES"
        ),
        "Angiogenesis"=c(
            "PID VEGF VEGFR PATHWAY"
        ),
        "Deregulating cellular energetics"=c(
            "REACTOME GLYCOLYSIS",
            "REACTOME RESPIRATORY ELECTRON TRANSPORT",
            "REACTOME PEROXISOMAL LIPID METABOLISM"
        ),
        "Avoiding immune destruction"=c(
            "BIOCARTA CTLA4 PATHWAY"
        ),
        "Tumor-promoting inflammation"=c(
            "BIOCARTA INFLAM PATHWAY"
        )
    )
