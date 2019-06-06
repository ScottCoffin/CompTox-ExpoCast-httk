#' Parameterize_PBTK
#' 
#' This function initializes the parameters needed in the functions solve_pbtk,
#' calc_css, and others using the multiple compartment model.
#' 
#' 
#' @param chem.name Either the chemical name or the CAS number must be
#' specified. 
#' @param chem.cas Either the chemical name or the CAS number must be
#' specified. 
#' @param species Species desired (either "Rat", "Rabbit", "Dog", "Mouse", or
#' default "Human").
#' @param default.to.human Substitutes missing animal values with human values
#' if true (hepatic intrinsic clearance or fraction of unbound plasma).
#' @param tissuelist Specifies compartment names and tissues groupings.
#' Remaining tissues in tissue.data are lumped in the rest of the body.
#' However, solve_pbtk only works with the default parameters.
#' @param force.human.clint.fup Forces use of human values for hepatic
#' intrinsic clearance and fraction of unbound plasma if true.
#' @param clint.pvalue.threshold Hepatic clearance for chemicals where the in
#' vitro clearance assay result has a p-values greater than the threshold are
#' set to zero.
#' @param adjusted.Funbound.plasma Returns adjusted Funbound.plasma when set to
#' TRUE along with parition coefficients calculated with this value.
#' @param regression Whether or not to use the regressions in calculating
#' partition coefficients.
#' @param suppress.messages Whether or not the output message is suppressed.
#' @return
#' @param minimum.Funbound.plasma Monte Carlo draws less than this value are set 
#' equal to this value (default is 0.0001 -- half the lowest measured Fup in our
#' dataset).
#' 
#' \item{BW}{Body Weight, kg.} \item{Clmetabolismc}{Hepatic Clearance, L/h/kg
#' BW.} \item{Fgutabs}{Fraction of the oral dose absorbed, i.e. the fraction of
#' the dose that enters the gutlumen.} \item{Funbound.plasma}{Fraction of
#' plasma that is not bound.} \item{Fhep.assay.correction}{The fraction of
#' chemical unbound in hepatocyte assay using the method of Kilford et al.
#' (2008)} \item{hematocrit}{Percent volume of red blood cells in the blood.}
#' \item{Kgut2pu}{Ratio of concentration of chemical in gut tissue to unbound
#' concentration in plasma.} \item{kgutabs}{Rate that chemical enters the gut
#' from gutlumen, 1/h.} \item{Kkidney2pu}{Ratio of concentration of chemical in
#' kidney tissue to unbound concentration in plasma.} \item{Kliver2pu}{Ratio of
#' concentration of chemical in liver tissue to unbound concentration in
#' plasma.} \item{Klung2pu}{Ratio of concentration of chemical in lung tissue
#' to unbound concentration in plasma.} \item{Krbc2pu}{Ratio of concentration
#' of chemical in red blood cells to unbound concentration in plasma.}
#' \item{Krest2pu}{Ratio of concentration of chemical in rest of body tissue to
#' unbound concentration in plasma.} \item{million.cells.per.gliver}{Millions
#' cells per gram of liver tissue.} \item{MW}{Molecular Weight, g/mol.}
#' \item{Qcardiacc}{Cardiac Output, L/h/kg BW^3/4.} \item{Qgfrc}{Glomerular
#' Filtration Rate, L/h/kg BW^3/4, volume of fluid filtered from kidney and
#' excreted.} \item{Qgutf}{Fraction of cardiac output flowing to the gut.}
#' \item{Qkidneyf}{Fraction of cardiac output flowing to the kidneys.}
#' \item{Qliverf}{Fraction of cardiac output flowing to the liver.}
#' \item{Rblood2plasma}{The ratio of the concentration of the chemical in the
#' blood to the concentration in the plasma from available_rblood2plasma.}
#' \item{Vartc}{Volume of the arteries per kg body weight, L/kg BW.}
#' \item{Vgutc}{Volume of the gut per kg body weight, L/kg BW.}
#' \item{Vkidneyc}{Volume of the kidneys per kg body weight, L/kg BW.}
#' \item{Vliverc}{Volume of the liver per kg body weight, L/kg BW.}
#' \item{Vlungc}{Volume of the lungs per kg body weight, L/kg BW.}
#' \item{Vrestc}{ Volume of the rest of the body per kg body weight, L/kg BW.}
#' \item{Vvenc}{Volume of the veins per kg body weight, L/kg BW.} 
#' @author John Wambaugh and Robert Pearce
#' @references Kilford, P. J., Gertz, M., Houston, J. B. and Galetin, A.
#' (2008). Hepatocellular binding of drugs: correction for unbound fraction in
#' hepatocyte incubations using microsomal binding or drug lipophilicity data.
#' Drug Metabolism and Disposition 36(7), 1194-7, 10.1124/dmd.108.020834.
#' @keywords Parameter
#' @examples
#' 
#' 
#'  parameters <- parameterize_pbtk(chem.cas='80-05-7')
#' 
#'  parameters <- parameterize_pbtk(chem.name='Bisphenol-A',species='Rat')
#' 
#'  # Change the tissue lumping (note, these model parameters will not work with our current solver):
#'  compartments <- list(liver=c("liver"),fast=c("heart","brain","muscle","kidney"),
#'                       lung=c("lung"),gut=c("gut"),slow=c("bone"))
#'  parameterize_pbtk(chem.name="Bisphenol a",species="Rat",default.to.human=TRUE,
#'                    tissuelist=compartments) 
#'  
#'  
#' 
#' @export parameterize_pbtk
parameterize_pbtk <- function(chem.cas=NULL,
                              chem.name=NULL,
                              species="Human",
                              default.to.human=F,
                              tissuelist=list(
                                liver=c("liver"),
                                kidney=c("kidney"),
                                lung=c("lung"),
                                gut=c("gut")),
                              force.human.clint.fup = F,
                              clint.pvalue.threshold=0.05,
                              adjusted.Funbound.plasma=T,
                              regression=T,
                              suppress.messages=F,
                              minimum.Funbound.plasma=0.0001,
                              Caco2.options = list(Caco2.Pab.default = 2,
                                                   Caco2.Fgut = TRUE,
                                                   Caco2.Fabs = TRUE)
                              )
{
  physiology.data <- physiology.data
# Look up the chemical name/CAS, depending on what was provide:
  out <- get_chem_id(chem.cas=chem.cas,chem.name=chem.name)
  chem.cas <- out$chem.cas
  chem.name <- out$chem.name
   
  if(class(tissuelist)!='list') stop("tissuelist must be a list of vectors.") 
  # Clint has units of uL/min/10^6 cells
  Clint.db <- try(get_invitroPK_param("Clint",species,chem.CAS=chem.cas),silent=T)
  # Check that the trend in the CLint assay was significant:
  Clint.pValue <- try(get_invitroPK_param("Clint.pValue",species,chem.CAS=chem.cas),silent=T)
  if ((class(Clint.db) == "try-error" & default.to.human) || force.human.clint.fup) 
  {
    Clint.db <- try(get_invitroPK_param("Clint","Human",chem.CAS=chem.cas),silent=T)
    Clint.pValue <- try(get_invitroPK_param("Clint.pValue","Human",chem.CAS=chem.cas),silent=T)
    warning(paste(species,"coerced to Human for metabolic clearance data."))
  }
  if (class(Clint.db) == "try-error") stop("Missing metabolic clearance data for given species. Set default.to.human to true to substitute human value.")
  # Check if clint is a point value or a distribution, if a distribution, use the median:
  if (nchar(Clint.db) - nchar(gsub(",","",Clint.db))==3) 
  {
    Clint.dist <- Clint.db
    Clint<- as.numeric(strsplit(Clint.db,",")[[1]][1])
    Clint.pValue <- as.numeric(strsplit(Clint.db,",")[[1]][4])
    if (!suppress.messages) warning("Clint is provided as a distribution.")
  } else {
    Clint <- Clint.db
    Clint.dist <- NA
  }
  if (!is.na(Clint.pValue) & Clint.pValue > clint.pvalue.threshold) Clint  <- 0

  
# Predict the PCs for all tissues in the tissue.data table:
  schmitt.params <- parameterize_schmitt(chem.cas=chem.cas,
                                         species=species,
                                         default.to.human=default.to.human,
                                         force.human.fup=force.human.clint.fup,
                                         suppress.messages=T,
                                         minimum.Funbound.plasma=minimum.Funbound.plasma)
  PCs <- predict_partitioning_schmitt(parameters=schmitt.params,
                                      species=species,
                                      adjusted.Funbound.plasma=adjusted.Funbound.plasma,
                                      regression=regression,
                                      minimum.Funbound.plasma=minimum.Funbound.plasma)
# Get_lumped_tissues returns a list with the lumped PCs, vols, and flows:
  lumped_params <- lump_tissues(PCs,tissuelist=tissuelist,species=species)
  
# Check to see if we should use the in vitro fup assay correction:  
  if (adjusted.Funbound.plasma)
  {
    fup <- schmitt.params$Funbound.plasma
    warning('Funbound.plasma adjusted for in vitro partioning (Pearce, 2017). Set adjusted.Funbound.plasma to FALSE to use original value.')
  } else fup <- schmitt.params$unadjusted.Funbound.plasma

# Restrict the value of fup:
  if (fup < minimum.Funbound.plasma) fup <- minimum.Funbound.plasma

  Fgutabs <- try(get_invitroPK_param("Fgutabs",species,chem.CAS=chem.cas),silent=T)
  if (class(Fgutabs) == "try-error") Fgutabs <- 1
    
  
 # Check the species argument for capitilization problems and whether or not it is in the table:  
  if (!(species %in% colnames(physiology.data)))
  {
    if (toupper(species) %in% toupper(colnames(physiology.data)))
    {
      phys.species <- colnames(physiology.data)[toupper(colnames(physiology.data))==toupper(species)]
    } else stop(paste("Physiological PK data for",species,"not found."))
  } else phys.species <- species

# Load the physiological parameters for this species
  this.phys.data <- physiology.data[,phys.species]
  names(this.phys.data) <- physiology.data[,1]
  
  MW <- get_physchem_param("MW",chem.CAS=chem.cas) #g/mol
  pKa_Donor <- suppressWarnings(get_physchem_param("pKa_Donor",chem.CAS=chem.cas)) # acid dissociation constants
  pKa_Accept <- suppressWarnings(get_physchem_param("pKa_Accept",chem.CAS=chem.cas)) # basic association cosntants
  Pow <- 10^get_physchem_param("logP",chem.CAS=chem.cas) # Octanol:water partition coeffiecient

  outlist <- list()
   # Begin flows:
  #mL/min/kgBW converted to L/h/kgBW:
  QGFRc <- this.phys.data["GFR"]/1000*60 
  Qcardiacc = this.phys.data["Cardiac Output"]/1000*60 
  flows <- unlist(lumped_params[substr(names(lumped_params),1,1) == 'Q'])

  outlist <- c(outlist,c(
    Qcardiacc = as.numeric(Qcardiacc),
    flows[!names(flows) %in% c('Qlungf','Qtotal.liverf')],
    Qliverf= flows[['Qtotal.liverf']] - flows[['Qgutf']],
    Qgfrc = as.numeric(QGFRc))) 
  # end flows  
  
  # Begin volumes
  # units should be L/kgBW  
  Vartc = this.phys.data["Plasma Volume"]/(1-this.phys.data["Hematocrit"])/2/1000 #L/kgBW
  Vvenc = this.phys.data["Plasma Volume"]/(1-this.phys.data["Hematocrit"])/2/1000 #L/kgBW

  outlist <- c(outlist,
    Vartc = as.numeric(Vartc),
    Vvenc = as.numeric(Vvenc),
    lumped_params[substr(names(lumped_params),1,1) == 'V'],
    lumped_params[substr(names(lumped_params),1,1) == 'K'])
  
  
# Create the list of parameters:
  BW <- this.phys.data["Average BW"]
  hematocrit = this.phys.data["Hematocrit"]
  outlist <- c(outlist,list(BW = as.numeric(BW),
    kgutabs = 2.18, # 1/h
    Funbound.plasma = fup, # unitless fraction
    Funbound.plasma.dist = schmitt.params$Funbound.plasma.dist,
    hematocrit = as.numeric(hematocrit), # unitless ratio
    MW = MW, #g/mol
    Pow = Pow,
    pKa_Donor=pKa_Donor,
    pKa_Accept=pKa_Accept,
    MA=schmitt.params[["MA"]]))
  
  # Correct for unbound fraction of chemical in the hepatocyte intrinsic clearance assay (Kilford et al., 2008)
 outlist <- c(outlist,list(
              Fhep.assay.correction=calc_fu_hep(schmitt.params$Pow,
                pKa_Donor=schmitt.params$pKa_Donor,
                pKa_Accept=schmitt.params$pKa_Accept)))  # fraction 

  outlist <- c(outlist,
    list(Clint=Clint,
         Clint.dist = Clint.dist,
         Clmetabolismc= as.numeric(calc_hepatic_clearance(hepatic.model="unscaled",
                          parameters=list(
                            Clint=Clint, #uL/min/10^6 cells
                            Funbound.plasma=fup, # unitless fraction
                            Fhep.assay.correction=outlist$Fhep.assay.correction, 
                            million.cells.per.gliver= 110, # 10^6 cells/g-liver
                            liver.density= 1.05, # g/mL
                            Dn=0.17,BW=BW,
                            Vliverc=lumped_params$Vliverc, #L/kg
                            Qtotal.liverc=(lumped_params$Qtotal.liverc)/1000*60),
                          suppress.messages=T)), #L/h/kg BW
         million.cells.per.gliver=110, # 10^6 cells/g-liver
         liver.density=1.05 # g/mL
         )) 
  
  # Calculate Fgutabs
  # Caco-2 Pab:
  Caco2.Pab.db <- try(get_invitroPK_param("Caco2.Pab", species = "Human", chem.CAS = chem.cas), silent = T)
  if (class(Caco2.Pab.db) == "try-error"){  
    Caco2.Pab.db <- Caco2.options$Caco2.Pab.default
    warning(paste0("Default value of ", Caco2.options$Caco2.Pab.default, " used for Caco2 permeability."))
  }
  # Check if Caco2 a point value or a distribution, if a distribution, use the median:
  if (nchar(Caco2.Pab.db) - nchar(gsub(",","",Caco2.Pab.db)) == 2) 
  {
    Caco2.Pab.dist <- Caco2.Pab.db
    Caco2.Pab.point <- as.numeric(strsplit(Caco2.Pab.db,",")[[1]][1])
    if (!suppress.messages) warning("Clint is provided as a distribution.")
  } else {
    Caco2.Pab.point <- as.numeric(Caco2.Pab.db)
    Caco2.Pab.dist <- NA
  }
  gut.params <- list("cl_us" = outlist$Clmetabolismc, "BW" = BW, "Caco2.Pab" = Caco2.Pab.point)
  if(Caco2.options$Caco2.Fgut == FALSE){
    fgut.oral <- 1
  }else{
    fgut.oral <- calc_fgut.oral(Params = gut.params, species = species)
  }
  if(Caco2.options$Caco2.Fabs == FALSE){
    fabs.oral <- try(get_invitroPK_param("Fgutabs",species,chem.CAS=chem.cas),silent=T)
    if (class(fabs.oral) == "try-error") fabs.oral <- 1
  }else{
    fabs.oral <- calc_fabs.oral(Params = gut.params, species = "Human") # only calculable for human, assume the same across species
  }
  Fgutabs <- fabs.oral * fgut.oral

  outlist[["Fgutabs"]] <- Fgutabs
  
  if (adjusted.Funbound.plasma) 
  {
    outlist["Funbound.plasma.adjustment"] <- schmitt.params$Funbound.plasma.adjustment
  } else outlist["Funbound.plasma.adjustment"] <- NA
   
    outlist <- c(outlist,
      Rblood2plasma=available_rblood2plasma(chem.cas=chem.cas,
        species=species,
        adjusted.Funbound.plasma=adjusted.Funbound.plasma))
        
  return(outlist[sort(names(outlist))])
}
