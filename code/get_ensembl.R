# Function get_ensembl pulls ensembl datasets for each clade
get_ensembl <- function(clade){
  filename <- paste("../data/",
                    if(clade == "metazoa"){
                      "metazoa_ensembl.RData"
                    } else if(clade == "plants"){
                      "plants_ensembl.RData"
                    } else if(clade == "protists"){
                      "protists_ensembl.RData"
                    } else if(clade == "fungi"){
                      "fungi_ensembl.RData"
                    } else{
                      "chordata.ensembl.RData"
                    }, sep = "")
  if(!file.exists(filename)){
    ensembl <- useEnsembl(
      biomart = 
        if(clade == "metazoa"){
          "metazoa_mart"
        } else if(clade == "plants"){
          "plants_mart"
        } else if(clade == "protists"){
          "protists_mart"
        } else if(clade == "fungi"){
          "fungi_mart"
        } else{
          "ensembl"
        },
      host = 
        if(clade == "metazoa"){
          "metazoa.ensembl.org"
        } else if(clade == "plants"){
          "plants.ensembl.org"
        } else if(clade == "protists"){
          "protists.ensembl.org"
        } else if(clade == "fungi"){
          "fungi.ensembl.org"
        } else{
          "ensembl.org"
        })
    saveRDS(ensembl, file = filename)
    return(ensembl)
  }
  readRDS(filename)
}