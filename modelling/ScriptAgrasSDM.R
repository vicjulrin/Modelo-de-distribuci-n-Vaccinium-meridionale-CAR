# Organizar entorno de trabajo ####
## Cargar librerias necesarias para el análisis ####
# Lista de librerias necesarias
packages_list<-list("magrittr", "dplyr", "plyr", "this.path", "ggplot2", "raster", "terra", "sf", 
                    "maps", "tools", "spThin", "ENMeval","ggdendro", "tibble", "data.table", "pbapply", "openxlsx",
                    "future", "future.apply", "progressr", "performance", "igraph", "tidyr", "ggrepel", "MASS",
                    "pdp", "dismo", "ppcor", "gtools", "ggpubr", "gridExtra", "rlang", "rJava", "rrapply", "snow", "gtools", "pROC"
)

## Revisar e instalar librerias necesarias
packagesPrev<- .packages(all.available = TRUE)
lapply(packages_list, function(x) {   if ( ! x %in% packagesPrev ) { install.packages(x, force=T)}    })

## Cargar librerias
lapply(packages_list, library, character.only = TRUE)

## Establecer directorio de trabajo ####
dir_work<- this.path::this.path() %>% dirname()



# Definicion de parametros generales ####

## Establecer area M ####
areaM<-  terra::rast(file.path(dir_work, "areaM_convexBiomasAndinos.tif")) %>% setNames("sp")



# ruta archivo area de estudio M
path_areaM<- file.path(dir_work, "areaM_convexBiomasAndinos.tif")

# ruta folder donde hay variables ambeintales
folder_vars_areaM<- "D:/Repositorios/biom_rincon/biomodelos-sdm/modelling/agras_example_test/areaM_buffer1grado/Set_1"

# ruta folder variables ambientales
path_env_folder<- file.path(dir_work, "env_vars2")


# variables que deben mantenerse en pruebas de multicolinealidad
vars_predilection<-  c("topographic_earthenv_Elevation", "worldclim_bio12_Annual_Precipitation", "CLC_CLC2_Arbustales_rep_apot0.5km") %>% 
  {data.frame(Var= ., sort_pred=  seq_along(.))}

# numero de datos aleatorios - pseudo ausencias
Max.Bg<- 10000

# umbral de correlaicon para pruebas de multicolinealidad 
cor_threshold<- 0.65 

# umbral de probabilidad modelo de distribucion
model_threshold<- 0.75


# parametros prueba modelos 
model_options <- data.frame(
  fc= c("L", "Q", "P", "T", "H", "LQ", "LH", "LP", "LT", "QH", "QP", "QT", "HP", "HT", "PT",
        "LQH", "LQP", "LQT", "LPT", "LHT", "QHP", "QHT", "QPT", "HPT",
        "LQHP", "LQHT", "LQPT", "LHPT", "QHPT", "LQHPT"),
  fc_k= c("l", "q", "p", "t", "h", "lq", "lh", "lp", "lt", "qh", "qp", "qt", "ph", "th", "pt",
          "lqh", "lqp", "lqt", "lpt", "lph", "qph", "qth", "qpt", "pth",
          "lqph", "lqth", "lqpt", "lpth", "qpth", "lqpth")
)

models_input <- expand.grid(c("L", "Q", "P", "T", "H",
                              "LQ", "LH", "LP", "LT", "QH", "QP", "QT","HP", "HT", "PT",
                              "LQH", "LQP", "LQT", "LPT", "LHT", "QHP", "QHT", "QPT", "HPT",
                              "LQHP", "LQHT", "LQPT", "LHPT", "QHPT", "LQHPT"),
                            seq(0.5, 6, 0.5)) %>% setNames(c("fc", "rm")) %>% 
  dplyr::mutate(model= paste0(rm, "_", tolower(fc))) %>% list(model_options) %>% plyr::join_all()


# folder donde guardar los resultados generales
folder_out<- "K:/Unidades compartidas/Solicitudes/Agras/7rev_CAR/bioma"; dir.create(folder_out, recursive = T)

# Nombre de folder donde guardar los resultados detallados
sp_name<- "agras_models"

# Cargar datos de aprovechamiento por veredas de interes
info_aprovechamiento <- read.xlsx(file.path(dir_work, "info_aprovechamiento.xlsx"))


# Definicion de variables - area de estudio ####

## Establecer covariables ####
dir_layers<- file.path(dir_work, "env_vars2") # folder donde se almacenan todas las covariables (formato tif o asc)
env_layers<- list.files( dir_layers, "\\.asc$", recursive = T, full.names = T)

### Cargar covariables
stac_layers<- terra::rast(env_layers) %>% terra::mask(areaM) 
names(stac_layers)<-  gsub(",", ".", names(stac_layers)); names(stac_layers)<-  gsub("-", "",  names(stac_layers))

# remover covariables que no dan informacion en el area de estudio
stac_names<- terra::freq(stac_layers) %>% dplyr::filter(!is.na(value)) %>% dplyr::filter(!value==0) %>% {names(stac_layers)[unique(.$layer)]}
envMstack<- stac_layers[[stac_names]]

# reescribir covariables para optimizar el analisis
dir.create(folder_vars_areaM, showWarnings = F, recursive = T)
files_envM<- list.files( folder_vars_areaM, "\\.asc$", recursive = T, full.names = F) %>% tools::file_path_sans_ext()

for (i in 1:nlyr(envMstack)) {  name_layer<- names(envMstack)[i] ;
if(!(name_layer %in% files_envM)){
  
  envMi <- envMstack[[i]]
  
    value <- values(envMi, na.rm = TRUE) %>% max()
    y <- ndec(x = value)
    
    if (y == 0) { datTyp <- "INT2S"; decinum <- 0 }
    if (y >= 1) { datTyp <- "FLT4S"; decinum <- 3 }
    
    writeRaster(x = if(is.factor(envMi)){envMi}else{round(envMi, digits = decinum)},
                filename = file.path(folder_vars_areaM, paste0( names(envMstack[[i]]), ".asc")), overwrite = T, NAflag = -9999, datatype = datTyp)
  }   }

env.Mfiles <- list.files(folder_vars_areaM, ".asc$", recursive = T, full.names = T)
env.M <- terra::rast(env.Mfiles)[[stac_names]]


# Organizacion de registros ####

## Cargar datos ####
raster_base<- terra::rast(areaM)
raster_ids<- raster_base %>%  terra::setValues(seq(ncell(.)))
data_agras<- read.csv(file.path(dir_work, "BDmodelo_CAR.csv"))

## Espacializacion de datos y eliminacion de duplicados ####
data_agras_centroids<- data_agras %>% 
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs= 4326) %>%
  {dplyr::mutate(., id= terra::extract(raster_ids, .)[,2])} %>% dplyr::select(c("id")) %>% 
  dplyr::filter(!duplicated(id)) %>% st_drop_geometry() %>% cbind(terra::xyFromCell(raster_base, .$id))

raster_occurrences<- raster_base; raster_occurrences[data_agras_centroids$id]<- 1

data_occs<- terra::cells(raster_occurrences) %>% {cbind(terra::xyFromCell(raster_occurrences, .), env.M[.])} %>% 
  as.data.frame() %>%  dplyr::rename("longitude"="x", "latitude"="y") %>% na.omit()

# generar envMstack ajustado - sin variables vacias
valid_vars <- sapply(data_occs, function(x) {if (is.factor(x)) {return(any(levels(x) == "1" & x == "1"))} else {T}}) %>% { Filter(function(x) isTRUE(x), .) } %>% names()
data_occs<- data_occs %>% dplyr::select(valid_vars)
envMstack_adjust <- envMstack[[  names(envMstack) %>% {.[. %in% valid_vars]}   ]]


## Generacion de datos de entrenamiento y evaluacion ####
sp_occs<- data_occs %>% st_as_sf(coords = c("longitude", "latitude")) %>% {terra::extract(raster_ids, .)[,2]}

# definir bg - background data - pseduo ausencias 
ids_background<- cells(areaM) %>% {.[!. %in% sp_occs]}
sample_background<- sample(x = ids_background, size = Max.Bg, replace = F)
data_Sbg<- env.M
Sbg <- cbind(terra::xyFromCell(env.M, sample_background), env.M[sample_background]) %>% as.data.frame() %>% 
  dplyr::rename("longitude"="x", "latitude"="y") %>% na.omit()


# Calibracion del modelo - Seleccion iterativa de variables  ####

# organizar data
vars_data_cor<- envMstack_adjust[terra::cells(areaM)]  %>% mutate_if(is.factor, ~ as.numeric(as.character(.)))
vars_test <- names(vars_data_cor)

## loop iterativo multicolinealidad - importancia de variables
list_hclust_model<-  list()    
vars_test_vif_model<- vars_test

for(j in seq_along(vars_test_vif_model)){ print(paste0("j_", j))
  
  ## Test de multicolinealidad ####
  for(v in seq_along(vars_test_vif_model)){ print(paste0("v_", v))
    
    vars_data_test_model<- dplyr::select(vars_data_cor, vars_test_vif_model)
    
    ## Prueba de correlacion ####
    cordata_modelR<- cor(vars_data_test_model, method = "spearman") %>% as.data.frame.matrix()  %>% set_rownames(names(.))
    NACol<- names(which(rowSums(is.na(cordata_modelR)) > (ncol(cordata_modelR)/2) ))
    cordata_model<- cordata_modelR %>% {.[!names(.) %in% NACol,]} %>% {.[,!colnames(.) %in% NACol]}; 
    cordata_model[is.na(cordata_model)]<-0
    
    corhclust_model <- hclust(as.dist(1-abs(cordata_model)))
    group_covars_model<- cutree(corhclust_model, h = 1-cor_threshold) %>% as.data.frame %>% rownames_to_column("Var") %>%
      data.table::setnames(".", "group")
    
    ## Estimacion VIF ####
    vif_data_model <- ginv(as.matrix(cordata_model)) %>% diag() %>% setNames(names(cordata_model)) %>% as.data.frame() %>% rownames_to_column("Var") %>% 
      setNames(c("Var", "VIF")) %>% list(group_covars_model,vars_predilection) %>% plyr::join_all() %>% arrange(group, sort_pred, VIF)
    
    list_hclust_model[[length(list_hclust_model)+1]]<- list(corhclust=corhclust_model, vif_data= vif_data_model, cor_data= cordata_model)
    
    check_duplicated<- any(duplicated(group_covars_model$group))
    
    if( check_duplicated ){
      
      vars_test_vif_model<- vif_data_model %>% dplyr::filter(!duplicated(group))  %>% {.$Var}
      
    } else {
      break
    }
  }
  
  
  
  
 ## Test de calibracion del modelo  ####
  tune_list<- list(fc = unique(models_input$fc), rm = as.numeric(unique(models_input$rm)))
  data_occs_test<- data_occs %>% dplyr::select(c("longitude", "latitude", vars_test_vif_model))  
  Sbg_test<- Sbg %>% dplyr::select(c("longitude", "latitude", vars_test_vif_model)) 
  
  proc_userval<- function(vars) {
    proc <- kuenm::kuenm_proc(vars$occs.val.pred, c(vars$bg.train.pred, vars$bg.val.pred))
    out <- data.frame(
      proc_auc_ratio = proc$pROC_summary[1],
      proc_pval = proc$pROC_summary[2], row.names = NULL
    )
    return(out)
  }
    
  eval_all <- ENMevaluate(
    occs = data_occs_test , bg = Sbg_test, partitions = c("block"),
    tune.args = tune_list, algorithm = "maxent.jar",
    doClamp = T, user.eval = proc_userval, parallel = T, numCores= 25
  )
  
  future::plan(sequential); gc()
  length_models<- length(eval_all@variable.importance)
  
  ### Importancia de variables entre modelos ####
  importance_data <- eval_all@variable.importance %>% plyr::rbind.fill() %>% dplyr::group_by(variable) %>% 
    dplyr::summarise(
      percent.contribution = sum(percent.contribution)/length_models,
      permutation.importance = sum(permutation.importance)/length_models,
    ) %>% arrange(-permutation.importance) %>% dplyr::rename(Var=variable)
  
  list_hclust_model[[length(list_hclust_model)]]$importance<- importance_data
  list_hclust_model[[length(list_hclust_model)]]$vif_data<- list(importance_data, vif_data_model) %>% plyr::join_all()
  
  if( !any(importance_data$permutation.importance<1)  ){
    break
  } else {
    ### Eliminacion de variables no informativas 1% ####
    vars_test_vif_model<- importance_data %>% as.data.frame() %>%  arrange(permutation.importance) %>% dplyr::filter(permutation.importance>=1) %>% {.$Var}
  }
}

# Generacion de modelos calibrados -  validación cruzada jackknife  ####

## Guardar/cargar modelo compilado de seleccion de variables ####
# saveRDS(eval_all, file.path(folder_out, paste0("eval_compiled", ".rds")))
# eval_all <- readRDS(file.path(folder_out, paste0("eval_compiled", ".rds")))

tune_list<- list(fc = unique(models_input$fc), rm = as.numeric(unique(models_input$rm)))
data_occs_test_jack<- data_occs %>% dplyr::select(c("longitude", "latitude", eval_all@variable.importance[[1]]$variable) ) 
Sbg_test_jack<- Sbg %>% dplyr::select(c("longitude", "latitude", eval_all@variable.importance[[1]]$variable)) 

eval_jack <- ENMevaluate(
  occs = data_occs_test_jack , bg = Sbg_test_jack,
  tune.args = tune_list, partitions = "jackknife", algorithm = "maxent.jar",
  doClamp = T, user.eval = proc_userval, parallel = T, numCores= 25
)

## Guardar/cargar modelo calibrado ####
# saveRDS(eval_jack, file.path(folder_out, paste0("eval_jack", ".rds")))
# eval_jack <- readRDS(file.path(folder_out, paste0("eval_jack", ".rds")))
# eval_jack<- eval_all



# Generar modelos - mapas proyectados ####
eval_results<- eval_jack@results %>% list(model_options) %>% plyr::join_all() %>% dplyr::mutate(model= paste0(fc, "_", gsub("\\.", "", rm)), fc= fc_k)  %>%  dplyr::relocate(c("model"), .before = 1)
eval_results<- eval_results

## Organizar archivos para ejecutar la proyeccion maxent ####
folder_results<- paste0(gsub("/", "_", sp_name)) # nombre folder de resultados
dir.create(folder_results)

### Folder que almacena los resultados de evaluacion enmeval ####
setwd(folder_out)
folder_eval_results_enmeval<- file.path(folder_out, folder_results, "/eval_results_enmeval");
dir.create(folder_eval_results_enmeval, showWarnings = T, recursive = T)
write.csv(eval_results, file.path(folder_eval_results_enmeval, "eval_models.csv"), row.names = F )
write.csv(eval_results, file.path(folder_eval_results_enmeval, "eval_models.csv"), row.names = F )
best_kuenm_style <- data.frame(Model = as.character(paste0("M_", eval_results$rm, "_F_", tolower(eval_results$fc), "_Set_1")), name_model= eval_results$model)
write.csv(best_kuenm_style, file.path(folder_eval_results_enmeval, "selected_models.csv"), row.names = F)

### Folder que almacena los datos de ocurrencia ####
folder_occurrences<- file.path(folder_out, folder_results, "occurrences");
dir.create(folder_occurrences, showWarnings = T, recursive = T)
occ_kuenm <- eval_jack@occs[, c("longitude", "latitude")] %>% dplyr::mutate(species= basename(folder_results) ) %>% dplyr::relocate(species, .before = 1)
write.csv(occ_kuenm, file.path(folder_occurrences, "occ_joint_kuenm.csv"), row.names = F)

### Folder que almacena las variables de proyeccion ####
folder_M_variables<- file.path(folder_out, folder_results, "M_variables/Set_1");
dir.create(folder_M_variables, showWarnings = F, recursive = T)

vars_eval<- eval_jack@variable.importance[[1]]$variable
lapply(vars_eval, function(z) {
  file.copy(file.path(folder_vars_areaM, paste0(z, ".asc")),
            file.path(folder_out, paste0(folder_results, "/M_variables/Set_1"), paste0(z, ".asc")), recursive = T)
})


### Folder donde se almacenaran los resultados de proyeccion ####
folder_results_predict <- file.path(folder_out, folder_results, "/final_models_enmeval");
dir.create(folder_results_predict, showWarnings = F)

## Estimar modelos de proyeccion maxent ####
setwd(folder_out)
kuenm::kuenm_mod(
  occ.joint = paste0(folder_results, "/occurrences/occ_joint_kuenm.csv"),
  M.var.dir = paste0(folder_results, "/M_variables"), 
  out.eval = paste0(folder_results, "/eval_results_enmeval"),
  batch = paste0(folder_results, "/final_models"), rep.n = 1, rep.type = "Bootstrap",
  jackknife = F, out.dir = paste0(folder_results, "/final_models_enmeval"),
  max.memory = 2000, out.format = "cloglog",
  project = F, 
  ext.type = "no_ext", 
  write.mess = FALSE,
  write.clamp = FALSE, 
  maxent.path = "D:/Repositorios/biom_rincon/biomodelos-sdm/modelling", 
  args = c(NULL, paste0("maximumbackground=", nrow(eval_jack@bg) )), wait = TRUE, run = TRUE
)


# Estimacion metricas modelos ajustados ###
best_kuenm_style <- data.frame(Model = as.character(paste0("M_", eval_results$rm, "_F_", tolower(eval_results$fc), "_Set_1")), name_model= eval_results$model)
info_models<- best_kuenm_style %>% dplyr::mutate(name_eval= names(eval_jack@models), folder= file.path(folder_out,sp_name, "final_models_enmeval", name_model) )


## Ajustar nombres de los modelos ####
for(i in seq(nrow(info_models))){
  old_name<- file.path(folder_out,sp_name, "final_models_enmeval", best_kuenm_style[i,"Model"])
  new_name<- file.path(folder_out,sp_name, "final_models_enmeval", best_kuenm_style[i,"name_model"])
  file.rename(old_name, new_name)
}





# Estimar metricas por modelo individual ####
list_models<- info_models %>% split(.$name_model)
list_models2<- list_models

options(future.globals.maxSize = 2000  * 1024^2)

plan(sequential)
plan(multisession, workers= 20)

with_progress(result <- {
  p <- progressor(along= seq_along(list_models2))
  
  result<-future_lapply( list_models2, function(x) {p();
    
 

### Cargar datos especificos del modelo ####
eval_model<- eval_jack@models[[x$name_eval]]
eval_importance_model<- eval_jack@variable.importance[[x$name_eval]]
eval_result_model<- eval_jack@results %>% dplyr::filter(tune.args %in% x$name_eval  ) %>% cbind(x) %>% dplyr::relocate(names(x), .before = 1)

data_train_model<- list(eval_model@presence, eval_model@absence) %>% plyr::rbind.fill() %>% mutate_if(is.integer, ~ as.factor(as.character(.))) %>% na.omit()

map_model <- terra::rast(file.path(x$folder, "agras_models.asc")) %>% setNames("map_model")
map_model_treshold<- map_model;map_model_treshold[map_model_treshold<model_threshold]<- NA



### Figura de inmportancia de variables - modelo ####
summ_importance_data_model<- eval_importance_model %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(ranking= mean(c(percent.contribution, permutation.importance))) %>% 
  arrange(ranking) %>% 
  as.data.frame() %>% dplyr::mutate(variable= factor(variable, levels = unique(.$variable)))

dataplot_impvars_model<-  summ_importance_data_model %>% 
  pivot_longer(cols = c(percent.contribution, permutation.importance), 
               names_to = "tipo", values_to = "aporte") %>% 
  dplyr::mutate(tipo=factor(tipo, levels= rev(c("percent.contribution", "permutation.importance"))))

plot_impvars_model<-     ggplot(dataplot_impvars_model, aes(x = variable)) +
  geom_bar(aes(y = aporte, fill = tipo), stat = "identity", alpha= 0.5)+
  scale_y_continuous( breaks = seq(0, 100, by = 10)) +
  coord_cartesian(ylim = c(0,100))+
  
  coord_flip() + 
  labs(x = "Variables", y = "Porcentaje", fill = "") +
  scale_fill_manual(values = setNames(c("#0000FE", "#00ADAC"), 
                                      c("percent.contribution", "permutation.importance")),
                    labels = setNames(c("Contribución", "Pérdida\npermutación"), 
                                      c("percent.contribution", "permutation.importance"))
  ) + theme_minimal()

## Analisis de dependencia parcial - atributos de variables #####

#### Depedencia parcial de combinaciones ####
vars_pdp_model<- summ_importance_data_model  %>% arrange(-permutation.importance)
data_combinations_model<- as.character(vars_pdp_model$variable) %>%  {combinations(n = length(.), r = 2, v = ., repeats.allowed = TRUE)} %>% as.data.frame() %>% 
  setNames(c("Var1", "Var2")  ) %>% dplyr::mutate(comb= paste0("comb", seq(nrow(.))))


list_combinations_model<- data_combinations_model %>% split(.$comb)

list_prev_model<- list()
list_pdp_vars_model <- pblapply( list_combinations_model, function(y) {
  
  vars_model<- sort(unique(as.character(unlist(y[,c("Var1", "Var2")]))))
  index_var<- paste0(vars_model, collapse = "_")
  
    if(index_var %in% names(list_prev_model) ){ 
      pdp_model<- list_prev_model[[index_var]]
    } else {  
      data_grid<- lapply(vars_model, function(x) {
        type_var<- class(data_train_model[,x])
        if(type_var=="numeric"){
          quantile(data_train_model[,x], probs = seq(0, 1, length.out = 20)) } else { unique(data_train_model[,x]) }
      }) %>%  {do.call(expand.grid, .)} %>% setNames(vars_model) %>% dplyr::distinct()
      
      pdp_model <- pdp::partial(eval_model, pred.var = vars_model, pred.fun = function(object, newdata) {predict(object, newdata)} , train = data_train_model, ice = F, pred.grid = data_grid,
                                prob = T,  progress = "text", plot= F, chull = F, approx = T)
      
      list_prev_model[[index_var]]<<- pdp_model
    }
    
    factor_columns <- names(pdp_model)[sapply(pdp_model, is.factor)]
    if(length(factor_columns)>0){
      pdp_model_factor<- pdp_model
      for( j in factor_columns){
        pdp_model_factor<- pdp_model_factor[pdp_model_factor[,j] =="1", ]
      }
    } else {pdp_model_factor<- pdp_model}
    
    pdp_list<- lapply(vars_model, function(y) {
      
      type_var<- class(pdp_model[,y])
      
      if(type_var=="factor"){
        
        other_factors<- factor_columns %>% {.[!. %in% y]}
        if(length(other_factors)>0){
          pdp_model_0 <- pdp_model
          for( o in other_factors){
            pdp_model_0<- pdp_model_0[pdp_model_0[,o] =="1", ]
          }
        } else {pdp_model_0<- pdp_model}
        
        pdp_var<- pdp_model_0 %>% dplyr::group_by(!!sym(y)) %>% dplyr::summarise(yhat= mean(yhat, na.rm=T))
        pdp_var<- pdp_var[pdp_var[,y]==1,]
        
      } else {
        pdp_var<- pdp_model_factor %>% dplyr::group_by(!!sym(y)) %>% dplyr::summarise(yhat= mean(yhat, na.rm=T))
      }
      
      as.data.frame(pdp_var)
      
    }) %>% setNames(vars_model)
    
    pdp_list
  
})


#### Resumir metricas de combinaciones pdp ####
list_pdp_vars_unlist_summ <- c(list(compiled_pdp=data_combinations_model), list_pdp_vars_model)






  
#### Plot dependencia parcial - modelo ####
list_plots_pdpvar_model<- pblapply(vars_pdp_model$variable, function(x){


data_pdp_var_model<- dplyr::filter(data_combinations_model, Var1 %in% x | Var2 %in% x) %>%   
  dplyr::mutate(Var_comb = if_else(Var1 == x, Var2, Var1)) %>% mutate(Var_comb = factor(Var_comb, levels = rev(levels(vars_pdp_model$variable)))) %>%
  arrange(Var_comb) %>% dplyr::mutate(comb= factor(comb, levels= unique(.$comb))) %>%  split(.$comb)


plots_pdp_var_model<- pblapply(data_pdp_var_model, function(y) {
  
  data_comb<- list_pdp_vars_model[[as.character(y$comb)]][[as.character(x)]] %>% setNames(c("var", "prob"))
  
  label_x<- if(y$Var1 == x){y$Var2}else{y$Var1}
  
  type_var<- class(data_comb[, "var"])
  
  plot_pdp_var<- {if(type_var == "factor"){

    ggplot(data_comb, aes(x = var, y = prob)) +
      geom_bar(stat = "identity", fill = "#77ACD1", alpha= 0.8) +
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.75, ymax = 1), fill = "yellow", alpha = 0.1)+
      theme_light()+xlab(label_x)+ylab(x)+ theme(text = element_text(size = 4))

    
  } else {
    ggplot()+
      geom_line(data = data_comb, aes(x = var, y = prob), color= "#77ACD1")+
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.75, ymax = 1), fill = "yellow", alpha = 0.1)+
      theme_light()+xlab(label_x)+ylab(x)+ theme(text = element_text(size = 4))
  }}
  plot_pdp_var
  
})

compiled_plot_var<- ggpubr::ggarrange(plotlist = plots_pdp_var_model)

})


## Test de correlacion variables - modelo ####
valid_vars_model<- valid_vars %>% {.[!. %in% c("longitude", "latitude" )]}
layers_model<- terra::rast(env.Mfiles[grepl(paste(valid_vars, collapse= "|"), env.Mfiles)]) %>% {c(map_model, .)}
stac_data_model<- layers_model %>% {.[cells(map_model)]}

list_predict_var_model<- pblapply(names(stac_data_model)[2:ncol(stac_data_model)], function(x) {
  test_var<- lm(as.formula( paste0("map_model", "~", x) ), data = stac_data_model)  
  predict_var<- predict(test_var, stac_data_model) %>% as.data.frame() %>% setNames(x)
}) %>%  {do.call(cbind, .)}

partial_corr_model <- cor(list_predict_var_model) %>% as.data.frame.matrix()
data_cors_model<- partial_corr_model %>% dplyr::select(vars_pdp_model$variable)
cor_vars_model <- rownames(data_cors_model)[apply(data_cors_model, 1, function(row) any(abs(row) > cor_threshold))]

cordataR_model<- partial_corr_model[ cor_vars_model, cor_vars_model]
NACol<- names(which(rowSums(is.na(cordataR_model)) > (ncol(cordataR_model)/2) ))
cordata_model<- cordataR_model %>% {.[!names(.) %in% NACol,]} %>% {.[,!colnames(.) %in% NACol]}; cordata_model[is.na(cordata_model)]<-0

corhclust_model <- hclust(as.dist(1-abs(cordata_model))) 
cordend_model<-as.dendrogram(corhclust_model)
cordend_model_data <- dendro_data(cordend_model)

group_covars_model<- cutree(corhclust_model, h = 1-cor_threshold) %>% as.data.frame %>% rownames_to_column("Var") %>%
  data.table::setnames(".", "group")


#### Dendograma de correlacion variables - modelo ####
col1<- "#EBEBEB"; col2<- "white";
var_table_model <- with(cordend_model_data$labels, data.frame(y_center = x, y_min= x-0.5, y_max=x+0.5, Var = as.character(label), height = 1)) %>% 
  list(group_covars_model) %>% plyr::join_all() %>%  arrange(- y_center, group) %>%
  dplyr::rowwise() %>% 
  mutate(colVar = if_else(Var %in% vars_pdp_model$variable, "red", "black"),
         Var = if_else(Var %in% vars_pdp_model$variable, paste("*", Var , sep = ""), Var)
         ) %>% as.data.frame() %>% 
  mutate(change_flag = if_else(group != lag(group, default = first(group)), 1, 0),
         intercalated = cumsum(change_flag) %% 2) %>% dplyr::mutate(col= if_else(intercalated==1, col1, col2))

segment_data_model <- with(segment(cordend_model_data), data.frame(x = y, y = x, xend = yend, yend = xend, cor= 1-yend))

ggdendroPlot_Cor_map_model <-   ggplot()+
  annotate("rect", xmin = -0.05, xmax = 1.04, fill = var_table_model$col,ymin = var_table_model$y_min , ymax = var_table_model$y_max, alpha = 0.9 )+
  geom_segment(data= segment_data_model, aes(x = 1-x, y = y, xend = 1-xend, yend = yend, label= cor), size= 0.3) +
  scale_y_continuous(breaks = var_table_model$y_center,  labels = var_table_model$Var)+
  coord_cartesian(expand = F)+
  labs(x= "Correlation", y= "Variables") +
  geom_vline(xintercept = cor_threshold, linetype = "dashed", col= "red") +
  theme(legend.position =  "bottom", legend.key.width = unit(50, 'pt'),
        plot.margin = margin(t = 0, r = 0,  b = 0,l = 0),
        panel.grid.major = element_line(color = "gray"),
        axis.ticks.length   = unit(0.3, "mm"),
        text = element_text(size = 10),
        axis.text.y = element_text(colour = var_table_model$colVar)
  )

### Analisis de aprovechamiento - modelo ####
list_veredas_model<- info_aprovechamiento %>% split(.$Value)

rep_veredas_model<- pblapply(list_veredas_model, function(x) {
  
  rast_vereda <- terra::rast(x$path_rast) %>% terra::resample(map_model_treshold)
  rast_vereda_model <- map_model_treshold %>% terra::mask(rast_vereda)
  
  area_vereda_model <- terra::cellSize(rast_vereda_model, mask=T, unit="km") %>% values() %>% na.omit() %>% sum()
  
  rep_ver_model<- x %>% dplyr::mutate(area_vereda_model=area_vereda_model) %>% dplyr::mutate(rep_model_vereda= (area_vereda_model/area_vereda))
  rep_ver_model
}) %>% plyr::rbind.fill() %>% dplyr::select(- c("X3", "Value", "vereda", "path_rast"))

plot_aprovechamiento_model  <- ggplot(rep_veredas_model, aes(x = rep_model_vereda, y = Indice_aprov)) +
  geom_point(aes(label= Veredas , color= Departamento         ), shape= 17)+
  geom_text_repel(aes(label= Veredas , color= Departamento         ), size= 3, show_guide  = FALSE) +
  labs(x= "% Rep Modelo.", y= "Uso")   +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +  # Línea horizontal en y = 0.5
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  annotate("text", x = Inf, y = Inf, label = "I", hjust = 1.1, vjust = 1.1) +
  annotate("text", x = -Inf, y = Inf, label = "II", hjust = -0.1, vjust = 1.1) +
  annotate("text", x = -Inf, y = -Inf, label = "III", hjust = -0.1, vjust = - 1)+
  annotate("text", x = Inf, y = -Inf, label = "IV", hjust = 1.1, vjust = -1) +
  coord_cartesian(xlim= c(0,1), ylim= c(0,1))+
  scale_y_continuous(expand = c(0,0))


### Cortar resultados - area de interes - modelo ####
carcun_rast<- terra::rast("K:/Unidades compartidas/Solicitudes/Agras/carcun_raster.tif") %>% terra::resample(map_model)
mapmodel_carcun<- terra::mask(map_model, carcun_rast)
map_model_treshold_carcun<-  terra::mask(map_model_treshold, carcun_rast )
  


### Exportar resultados ####
folder_model<- x$folder
name_model<- gsub("\\.", "", x$name_model) 

# metricas modelo
saveRDS(eval_model, file.path(folder_model, paste0("ENMevaluation", "_", name_model, ".rds")))
openxlsx::write.xlsx(eval_result_model, file.path(folder_model, paste0("results_model_",name_model, ".xlsx")))
openxlsx::write.xlsx(eval_importance_model, file.path(folder_model, paste0("importance_data_",name_model, ".xlsx")))

# plots importancia de variables
ggsave(file.path(folder_model, paste0(  "impVar", name_model, ".png")), plot_impvars_model)
ggsave(file.path(folder_model, paste0(  "impVar", name_model, ".jpg")), plot_impvars_model)

# data correlaciones modelo final
openxlsx::write.xlsx(cordata_model, file.path(folder_model, paste0("impVar_cordata_",name_model, ".xlsx")))

# plot correlaciones modelo final
ggsave(file.path(folder_model, paste0(  "dendroCorVarsModel", name_model, ".png")), ggdendroPlot_Cor_map_model)
ggsave(file.path(folder_model, paste0(  "dendroCorVarsModel", name_model, ".jpg")), ggdendroPlot_Cor_map_model)

# plot analisis de aprovechamiento
openxlsx::write.xlsx(rep_veredas_model, file.path(folder_model, paste0("rep_veredas_model_",name_model, ".xlsx")))
ggsave(file.path(folder_model, paste0( "plot_uso_rep", name_model, ".png")), plot_aprovechamiento_model)
ggsave(file.path(folder_model, paste0("plot_uso_rep", name_model, ".jpg")), plot_aprovechamiento_model)

# raster resultados
terra::writeRaster(map_model, file.path(folder_model, paste0(  "full_", name_model, ".tif")),overwrite = T )
terra::writeRaster(map_model_treshold, file.path(folder_model, paste0(  "thresh_", gsub("\\.", "", model_threshold), "_" , name_model, ".tif") ),overwrite = T )

terra::writeRaster(mapmodel_carcun, file.path(folder_model, paste0(  "carcun_full_", name_model, ".tif")),overwrite = T )
terra::writeRaster(map_model_treshold_carcun, file.path(folder_model, paste0(  "carcun_thresh_", name_model, ".tif")),overwrite = T )


# metricas de dependencia parcial
openxlsx::write.xlsx(list_pdp_vars_unlist_summ, file.path(folder_model, paste0("impVar_pdp_attributes_",name_model, ".xlsx")))

# plot de dependencia parcial
pdf(  file.path(folder_model, paste0(  "impVar_pdp_attributes_plot_", name_model, ".pdf"))  )
for(plot in list_plots_pdpvar_model) { print(plot) }
dev.off(); dev.off()



rm(list = ls())
gc()




  }, future.chunk.size=1 )
  
  
})

# Estimar modelo ensamblado ####
## Cargar mapa de expertos ####
mapaExpertos<- terra::rast("K:/Unidades compartidas/Solicitudes/Agras/restricciones/mapaExpertosRestricciones.tif") %>% setNames("mapaExperto")

## Listar direcciones de modelos ####
ready_models<- list.files(folder_out, pattern = "impVar_pdp_attributes.*\\.xlsx$", recursive = T, full.names = T) 
folder_models <- dirname(ready_models) %>% unique()

folders_models_path <- pblapply(folder_models, function(x) {
  name_layer<- gsub("\\.", "", basename(x))
  list(folder= x, name_layer= name_layer) }) %>% setNames(sapply(., function(x) x$name_layer))





## Cargar stac de modelos ####
list_models_compiled <- pbsapply(folders_models_path, function(e) {
  file.path(e$folder, paste0(  "full_", e$name_layer, ".tif")) %>% terra::rast() %>% setNames(e$name_layer)
})   %>% terra::rast()


# organizar datos segun umbral
data_models_compiled<- c(mapaExpertos, list_models_compiled, raster_occurrences) %>% {.[cells(mapaExpertos)]} 
data_models_compiled_weights<- data_models_compiled %>% dplyr::mutate(mapaExperto= if_else(mapaExperto<0, mapaExperto/100000, 1))


eval_results2<- dplyr::filter(eval_results, !is.na(delta.AICc)) %>% arrange(delta.AICc)

list_evals<- eval_results2 %>% split(.$model)

score_expertos <-pblapply( list_evals , function(x) {
  
  mean_model_test<- data_models_compiled_weights %>% dplyr::select(c("mapaExperto", x$model)) 
  mean_model_tresh_test<- mean_model_test[mean_model_test[,x$model]>=0.75,]
  scores_expertos<- table(mean_model_tresh_test$mapaExperto) %>% {sapply(names(.), function(i) {   as.numeric(i)*(.[i])  }  )}
  score_ex<- scores_expertos[scores_expertos<0] %>% abs() %>% mean() %>% {if(is.nan(.)){0}else{.}}
  score_in<- scores_expertos[scores_expertos>0] %>% abs() %>% mean() %>% {if(is.nan(.)){0}else{.}}
  
  dplyr::mutate(x, score_ex= score_ex, score_in=score_in, score_expert= ( (x$AICc * ({if(score_ex==0){1}else{score_ex}}))/score_in ) )
} ) %>% plyr::rbind.fill() %>% arrange(score_expert)




## Estimacion iterativa mejor modelo ensamblado ####
test_ensamble <-pblapply( seq_along(score_expertos$model) , function(i) {print(i)
  
  better_models_test<- score_expertos[rev(i:1),]
  
  mean_model_tresh_test<- data_models_compiled_weights %>% dplyr::mutate(mean= rowMeans(dplyr::select(., better_models_test$model))) %>% 
    dplyr::select(c("mapaExperto", "mean")) %>% dplyr::filter(mean>=0.75)
  
  scores_expertos<- table(mean_model_tresh_test$mapaExperto) %>% {sapply(names(.), function(i) {   as.numeric(i)*(.[i])  }  )}
  score_ex<- scores_expertos[scores_expertos<0] %>% abs() %>% mean() %>% {if(is.nan(.)){0}else{.}}
  score_in<- scores_expertos[scores_expertos>0] %>% abs() %>% mean() %>% {if(is.nan(.)){0}else{.}}
  
  list(data_score= data.frame(index= i, score_ex=score_ex, score_in= score_in, score_expert= ( (mean(better_models_test$AICc) * ({if(score_ex==0){1}else{score_ex}}))/score_in ) )  , 
       better_models_test= better_models_test$model)
} ) 









results_test_ensamble <- purrr::map(test_ensamble, "data_score") %>% plyr::rbind.fill() %>% arrange( score_expert )
results_test_ensamble$ensamble_id<- seq(nrow(results_test_ensamble))



write.xlsx(score_expertos, "K:/Unidades compartidas/Solicitudes/Agras/7rev_CAR/bioma/agras_models/score_expertos.xlsx")



write.xlsx(results_test_ensamble, "K:/Unidades compartidas/Solicitudes/Agras/7rev_CAR/bioma/agras_models/ensamble_scores.xlsx")





i<- 2




library(future)
library(future.apply)
library(progressr)

plan(sequential)
plan(multisession, workers= 25)


with_progress(result <- {
  p <- progressor(along= seq(nrow(results_test_ensamble)) )
  
  result<-future_lapply( 3:nrow(results_test_ensamble), function(i) {
    




### Exportar resultados ####
folder_ensamble<- file.path(folder_out, "agras_models", "ensamble_models",  paste0("ensamble_model_", i) ); dir.create(folder_ensamble, recursive = T)
name_ensamble<- gsub("\\.", "", paste0("ensamble_model_", i) ) 


better_models_index<- test_ensamble[[results_test_ensamble[i,]$index]]
better_models_ensamble<- better_models_index$better_models_test

## Estadisticas mejor modelo ensamblado  ####

# Generacion modelo ensamblado
list_models_compiled <- pbsapply(folders_models_path, function(e) {
  file.path(e$folder, paste0(  "full_", e$name_layer, ".tif")) %>% terra::rast() %>% setNames(e$name_layer)
})   %>% terra::rast()

map_ensamble <- list_models_compiled[[ better_models_ensamble ]] %>% mean() %>% setNames("map_ensamble") 
map_ensamble_treshold<- map_ensamble;map_ensamble_treshold[map_ensamble_treshold<model_threshold]<- NA


# Compilar estadisticas del modelo ensamblado ####
folders_better_models_path<- folders_models_path[better_models_ensamble]

# resultados medios del modelo ensamblado
eval_results_models_ensamble<- score_expertos %>% dplyr::filter(model %in% better_models_ensamble)
eval_result_ensamble<- eval_results_models_ensamble %>%  select_if(is.numeric) %>% dplyr::select(-names(.)[grepl("score", names(.))]) %>% colMeans(na.rm=T) %>% {as.data.frame(t(.))} %>% 
  cbind(dplyr::select(better_models_index$data_score, -"index"))

# Importancia media del modelo ensamblado
eval_importance_ensamble<- pblapply(folders_better_models_path, function(x) {
  file.path(x$folder, paste0("importance_data_",x$name_layer, ".xlsx")) %>% read.xlsx() 
}) %>% plyr::rbind.fill() %>% dplyr::group_by(variable) %>% 
  dplyr::summarise(
    percent.contribution = sum(percent.contribution)/ length(folders_better_models_path),
    permutation.importance = sum(permutation.importance)/length(folders_better_models_path),
  ) %>% arrange(-permutation.importance) 



### Figura de inmportancia de variables - modelo ensamblado ####
summ_importance_data_ensamble<- eval_importance_ensamble %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(ranking= mean(c(percent.contribution, permutation.importance))) %>% 
  arrange(ranking) %>% 
  as.data.frame() %>% dplyr::mutate(variable= factor(variable, levels = unique(.$variable)))

dataplot_impvars_ensamble<-  summ_importance_data_ensamble %>% 
  pivot_longer(cols = c(percent.contribution, permutation.importance), 
               names_to = "tipo", values_to = "aporte") %>% 
  dplyr::mutate(tipo=factor(tipo, levels= rev(c("percent.contribution", "permutation.importance"))))

plot_impvars_ensamble<-     ggplot(dataplot_impvars_ensamble, aes(x = variable)) +
  geom_bar(aes(y = aporte, fill = tipo), stat = "identity", alpha= 0.5)+
  scale_y_continuous( breaks = seq(0, 100, by = 10)) +
  coord_cartesian(ylim = c(0,100))+
  
  coord_flip() + 
  labs(x = "Variables", y = "Porcentaje", fill = "") +
  scale_fill_manual(values = setNames(c("#0000FE", "#00ADAC"), 
                                      c("percent.contribution", "permutation.importance")),
                    labels = setNames(c("Contribución", "Pérdida\npermutación"), 
                                      c("percent.contribution", "permutation.importance"))
  ) + theme_minimal()


### Analisis de dependencia parcial - atributos de variables - modelo ensamblado #####
data_combinations_ensamble<- pblapply(folders_better_models_path, function(x) {
  path_pdp<- file.path(x$folder, paste0("impVar_pdp_attributes_",x$name_layer, ".xlsx"))
  model_pdp<- path_pdp %>% read.xlsx() %>% 
    dplyr::filter(Var1 %in% summ_importance_data_ensamble$variable | Var2 %in% summ_importance_data_ensamble$variable) %>% dplyr::distinct() %>% 
    dplyr::mutate(path_pdp= path_pdp) 
}) %>% plyr::rbind.fill() %>% dplyr::mutate(id_comb = paste0(Var1, "_", Var2)  ) %>% dplyr::mutate(id_comb = factor(id_comb, levels= unique(.$id_comb))  ) %>% 
  dplyr::mutate(id_comb= paste0("comb", as.numeric(id_comb)))



pdp_better_models_list<- data_combinations_ensamble %>% split(.$Var1)

list_pdp_vars_ensamble<- pblapply( pdp_better_models_list , function(y) { 
  
  list_Var2<- split(y, y$Var2)
  
  list_pdp_Var2<- pblapply(list_Var2, function(z) {
    
    name_vars<- dplyr::select(z, c("Var1", "Var2")) %>% dplyr::distinct() %>% unlist() %>% unique()
    
    data_pdp_list<- z %>% split(.$path_pdp) %>% {if(length(.)>1){.[1:2]}else {.}} %>%   lapply(function(z) {read.xlsx(z$path_pdp, z$comb, colNames = F, rowNames = F) }) %>% plyr::rbind.fill() 
    
    pdp_var2 <- lapply(seq(1, ncol(data_pdp_list), by = 2), function(i) { data_pdp_list[, i:min(i + 1, ncol(data_pdp_list))]
    }) %>% setNames(name_vars) %>%
      {lapply(names(.), function(h)  .[[h]] %>% setNames(c(h, "yhat")) %>% dplyr::group_by(!!sym(h)) %>% dplyr::summarise(yhat= mean(yhat, na.rm=T)) )} %>% setNames(name_vars)
    
  })
  list_pdp_Var2
})



list_pdp_vars_unlist_summ_ensamble<- c(list(compiled_pdp=data_combinations_ensamble %>%   dplyr::select(c("Var1", "Var2", "id_comb")) %>% dplyr::distinct() ), 
                                       unlist(list_pdp_vars_ensamble, recursive = F) %>% setNames( unique(data_combinations_ensamble$id_comb) ) )


#### Plot dependencia parcial - modelo ensamblado ####
vars_pdp_ensamble<- summ_importance_data_ensamble  %>% arrange(-permutation.importance)



list_plots_pdpvar_ensamble<- pblapply(vars_pdp_ensamble$variable, function(x){
  
  
  data_pdp_var_ensamble<- dplyr::filter(data_combinations_ensamble, Var1 %in% x | Var2 %in% x) %>%   dplyr::select(c("Var1", "Var2", "id_comb")) %>% dplyr::distinct() %>% 
    dplyr::mutate(Var_comb = if_else(Var1 == x, Var2, Var1)) %>% mutate(Var_comb = factor(Var_comb, levels = rev(levels(vars_pdp_ensamble$variable)))) %>%
    arrange(Var_comb) %>%  dplyr::mutate(id_comb= factor(id_comb, levels= unique(.$id_comb))) %>% split(.$id_comb)
  
  
  plots_pdp_var_ensamble<- pblapply(seq_along(data_pdp_var_ensamble), function(y) { print(y); y<- data_pdp_var_ensamble[[y]]
  
  data_comb<- list_pdp_vars_ensamble[[as.character(y$Var1)]][[y$Var2]][[as.character(x)]] %>% as.data.frame() %>% setNames(c("var", "prob"))
  
  label_x<- if(y$Var1 == x){y$Var2}else{y$Var1}
  
  type_var<- class(data_comb[, "var"])
  
  plot_pdp_var<- {if(type_var == "numeric"){
    
    ggplot()+
      geom_line(data = data_comb, aes(x = var, y = prob), color= "#77ACD1")+
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.75, ymax = 1), fill = "yellow", alpha = 0.1)+
      theme_light()+xlab(label_x)+ggtitle(x)+ theme(text = element_text(size = 4))+ylab("")

  } else {

    
    ggplot(data_comb, aes(x = var, y = prob)) +
      geom_bar(stat = "identity", fill = "#77ACD1", alpha= 0.8) +
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.75, ymax = 1), fill = "yellow", alpha = 0.1)+
      theme_light()+xlab(label_x)+ggtitle(x)+ theme(text = element_text(size = 4))+ylab("")
    
    
  }}
  plot_pdp_var
  
  })
  

})



col1<- list(list_plots_pdpvar_ensamble[[1]][[1]]+xlab(""), list_plots_pdpvar_ensamble[[2]][[1]]+xlab(""), list_plots_pdpvar_ensamble[[3]][[1]]+xlab(""), 
            list_plots_pdpvar_ensamble[[8]][[1]]+xlab(""),
            list_plots_pdpvar_ensamble[[5]][[1]])

col2<- list(list_plots_pdpvar_ensamble[[1]][[2]]+xlab("")+ggtitle(""), list_plots_pdpvar_ensamble[[2]][[2]]+xlab("")+ggtitle(""), list_plots_pdpvar_ensamble[[3]][[2]]+xlab("")+ggtitle(""), 
            list_plots_pdpvar_ensamble[[8]][[2]]+xlab("")+ggtitle(""),
            list_plots_pdpvar_ensamble[[5]][[2]]+ggtitle(""))

col3<- list(list_plots_pdpvar_ensamble[[1]][[3]]+xlab("")+ggtitle(""), list_plots_pdpvar_ensamble[[2]][[3]]+xlab("")+ggtitle(""), list_plots_pdpvar_ensamble[[3]][[3]]+xlab("")+ggtitle(""), 
            list_plots_pdpvar_ensamble[[8]][[3]]+xlab("")+ggtitle(""),
            list_plots_pdpvar_ensamble[[5]][[3]]+ggtitle(""))

col4<- list(list_plots_pdpvar_ensamble[[1]][[4]]+xlab("")+ggtitle(""), list_plots_pdpvar_ensamble[[2]][[4]]+xlab("")+ggtitle(""), list_plots_pdpvar_ensamble[[3]][[4]]+xlab("")+ggtitle(""), 
            list_plots_pdpvar_ensamble[[8]][[4]]+xlab("")+ggtitle(""),
            list_plots_pdpvar_ensamble[[5]][[4]]+ggtitle(""))






better_pdp<- ggarrange(plotlist = 
list(ggarrange(plotlist = col1, ncol = 1),
ggarrange(plotlist = col2, ncol = 1),
ggarrange(plotlist = col3, ncol = 1),
ggarrange(plotlist = col4, ncol = 1)),
ncol = 4
)

annotated_plot <- annotate_figure(better_pdp, 
                                  left = text_grob("Probabilidad predicha", rot = 90,  size = 6))


ggsave(file.path(folder_ensamble, paste0(  "better_pdp", name_ensamble, ".png")), annotated_plot)
ggsave(file.path(folder_ensamble, paste0(  "better_pdp", name_ensamble, ".jpg")), annotated_plot)

pdf(  file.path(folder_ensamble, paste0(  "better_pdp_attributes_plot_", name_ensamble, ".pdf")), width = 5, height = 5.5  )
print(annotated_plot)
dev.off(); dev.off()




## Test de correlacion variables - ensamblado ensamblado ####
valid_vars_ensamble<- valid_vars %>% {.[!. %in% c("longitude", "latitude" )]}
layers_ensamble<- terra::rast(env.Mfiles[grepl(paste(valid_vars, collapse= "|"), env.Mfiles)]) %>% {c(map_ensamble, .)}
stac_data_ensamble<- layers_ensamble %>% {.[cells(map_ensamble)]}

list_predict_var_ensamble<- pblapply(names(stac_data_ensamble)[2:ncol(stac_data_ensamble)], function(x) {
  test_var<- lm(as.formula( paste0("map_ensamble", "~", x) ), data = stac_data_ensamble)  
  predict_var<- predict(test_var, stac_data_ensamble) %>% as.data.frame() %>% setNames(x)
}) %>%  {do.call(cbind, .)}

partial_corr_ensamble <- cor(list_predict_var_ensamble) %>% as.data.frame.matrix()
data_cors_ensamble<- partial_corr_ensamble %>% dplyr::select(vars_pdp_ensamble$variable)
cor_vars_ensamble <- rownames(data_cors_ensamble)[apply(data_cors_ensamble, 1, function(row) any(abs(row) > cor_threshold))]

cordataR_ensamble<- partial_corr_ensamble[ cor_vars_ensamble, cor_vars_ensamble]
NACol<- names(which(rowSums(is.na(cordataR_ensamble)) > (ncol(cordataR_ensamble)/2) ))
cordata_ensamble<- cordataR_ensamble %>% {.[!names(.) %in% NACol,]} %>% {.[,!colnames(.) %in% NACol]}; cordata_ensamble[is.na(cordata_ensamble)]<-0

corhclust_ensamble <- hclust(as.dist(1-abs(cordata_ensamble))) 
cordend_ensamble<-as.dendrogram(corhclust_ensamble)
cordend_ensamble_data <- dendro_data(cordend_ensamble)

group_covars_ensamble<- cutree(corhclust_ensamble, h = 1-cor_threshold) %>% as.data.frame %>% rownames_to_column("Var") %>%
  data.table::setnames(".", "group")


#### Dendograma de correlacion variables - ensamblado ####
col1<- "#EBEBEB"; col2<- "white";
var_table_ensamble <- with(cordend_ensamble_data$labels, data.frame(y_center = x, y_min= x-0.5, y_max=x+0.5, Var = as.character(label), height = 1)) %>% 
  list(group_covars_ensamble) %>% plyr::join_all() %>%  arrange(- y_center, group) %>%
  dplyr::rowwise() %>% 
  mutate(colVar = if_else(Var %in% vars_pdp_ensamble$variable, "red", "black"),
         Var = if_else(Var %in% vars_pdp_ensamble$variable, paste("*", Var , sep = ""), Var)
  ) %>% as.data.frame() %>% 
  mutate(change_flag = if_else(group != lag(group, default = first(group)), 1, 0),
         intercalated = cumsum(change_flag) %% 2) %>% dplyr::mutate(col= if_else(intercalated==1, col1, col2))

segment_data_ensamble <- with(segment(cordend_ensamble_data), data.frame(x = y, y = x, xend = yend, yend = xend, cor= 1-yend))

ggdendroPlot_Cor_map_ensamble <-   ggplot()+
  annotate("rect", xmin = -0.05, xmax = 1.04, fill = var_table_ensamble$col,ymin = var_table_ensamble$y_min , ymax = var_table_ensamble$y_max, alpha = 0.9 )+
  geom_segment(data= segment_data_ensamble, aes(x = 1-x, y = y, xend = 1-xend, yend = yend, label= cor), size= 0.3) +
  scale_y_continuous(breaks = var_table_ensamble$y_center,  labels = var_table_ensamble$Var)+
  coord_cartesian(expand = F)+
  labs(x= "Correlation", y= "Variables") +
  geom_vline(xintercept = cor_threshold, linetype = "dashed", col= "red") +
  theme(legend.position =  "bottom", legend.key.width = unit(50, 'pt'),
        plot.margin = margin(t = 0, r = 0,  b = 0,l = 0),
        panel.grid.major = element_line(color = "gray"),
        axis.ticks.length   = unit(0.3, "mm"),
        text = element_text(size = 10),
        axis.text.y = element_text(colour = var_table_ensamble$colVar)
  )

### Analisis de aprovechamiento - ensamblado ####
list_veredas_ensamble<- info_aprovechamiento %>% split(.$Value)

rep_veredas_ensamble<- pblapply(list_veredas_ensamble, function(x) {
  
  rast_vereda <- terra::rast(x$path_rast) %>% terra::resample(map_ensamble_treshold)
  rast_vereda_ensamble <- map_ensamble_treshold %>% terra::mask(rast_vereda)
  
  area_vereda_ensamble <- terra::cellSize(rast_vereda_ensamble, mask=T, unit="km") %>% values() %>% na.omit() %>% sum()
  
  rep_ver_ensamble<- x %>% dplyr::mutate(area_vereda_ensamble=area_vereda_ensamble) %>% dplyr::mutate(rep_ensamble_vereda= (area_vereda_ensamble/area_vereda))
  rep_ver_ensamble
}) %>% plyr::rbind.fill() %>% dplyr::select(- c("X3", "Value", "vereda", "path_rast"))

plot_aprovechamiento_ensamble  <- ggplot(rep_veredas_ensamble, aes(x = rep_ensamble_vereda, y = Indice_aprov)) +
  geom_point(aes(label= Veredas , color= Departamento         ), shape= 17)+
  geom_text_repel(aes(label= Veredas , color= Departamento         ), size= 3, show_guide  = FALSE) +
  labs(x= "% Rep ensamblado.", y= "Uso")   +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +  # Línea horizontal en y = 0.5
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  annotate("text", x = Inf, y = Inf, label = "I", hjust = 1.1, vjust = 1.1) +
  annotate("text", x = -Inf, y = Inf, label = "II", hjust = -0.1, vjust = 1.1) +
  annotate("text", x = -Inf, y = -Inf, label = "III", hjust = -0.1, vjust = - 1)+
  annotate("text", x = Inf, y = -Inf, label = "IV", hjust = 1.1, vjust = -1) +
  coord_cartesian(xlim= c(0,1), ylim= c(0,1))+
  scale_y_continuous(expand = c(0,0))


### Cortar resultados - area de interes - ensamblado ####
carcun_rast<- terra::rast("K:/Unidades compartidas/Solicitudes/Agras/carcun_raster.tif") %>% terra::resample(map_ensamble)
map_ensamble_carcun<- terra::mask(map_ensamble, carcun_rast)
map_ensamble_treshold_carcun<-  terra::mask(map_ensamble_treshold, carcun_rast )


# metricas ensamblado
openxlsx::write.xlsx(eval_result_ensamble, file.path(folder_ensamble, paste0("results_ensamble_",name_ensamble, ".xlsx")))
openxlsx::write.xlsx(eval_importance_ensamble, file.path(folder_ensamble, paste0("importance_data_",name_ensamble, ".xlsx")))

# plots importancia de variables
ggsave(file.path(folder_ensamble, paste0(  "impVar", name_ensamble, ".png")), plot_impvars_ensamble)
ggsave(file.path(folder_ensamble, paste0(  "impVar", name_ensamble, ".jpg")), plot_impvars_ensamble)

# data correlaciones ensamblado final
openxlsx::write.xlsx(cordata_ensamble, file.path(folder_ensamble, paste0("impVar_cordata_",name_ensamble, ".xlsx")))

# plot correlaciones ensamblado final
ggsave(file.path(folder_ensamble, paste0(  "dendroCorVarsModel", name_ensamble, ".png")), ggdendroPlot_Cor_map_ensamble)
ggsave(file.path(folder_ensamble, paste0(  "dendroCorVarsModel", name_ensamble, ".jpg")), ggdendroPlot_Cor_map_ensamble)

# plot analisis de aprovechamiento
openxlsx::write.xlsx(rep_veredas_ensamble, file.path(folder_ensamble, paste0("rep_veredas_ensamble_",name_ensamble, ".xlsx")))
ggsave(file.path(folder_ensamble, paste0( "plot_uso_rep", name_ensamble, ".png")), plot_aprovechamiento_ensamble)
ggsave(file.path(folder_ensamble, paste0("plot_uso_rep", name_ensamble, ".jpg")), plot_aprovechamiento_ensamble)

# raster resultados
terra::writeRaster(map_ensamble, file.path(folder_ensamble, paste0(  "full_", name_ensamble, ".tif")),overwrite = T )
terra::writeRaster(map_ensamble_treshold, file.path(folder_ensamble, paste0(  "thresh_", gsub("\\.", "", model_threshold), "_" , name_ensamble, ".tif") ),overwrite = T )

terra::writeRaster(map_ensamble_carcun, file.path(folder_ensamble, paste0(  "carcun_full_", name_ensamble, ".tif")),overwrite = T )
terra::writeRaster(map_ensamble_treshold_carcun, file.path(folder_ensamble, paste0(  "carcun_thresh_", name_ensamble, ".tif")),overwrite = T )


# metricas de dependencia parcial
openxlsx::write.xlsx(list_pdp_vars_unlist_summ_ensamble, file.path(folder_ensamble, paste0("impVar_pdp_attributes_",name_ensamble, ".xlsx")))

# plot de dependencia parcial
pdf(  file.path(folder_ensamble, paste0(  "impVar_pdp_attributes_plot_", name_ensamble, ".pdf"))  )
for(plot in list_plots_pdpvar_ensamble) { print(plot) }
dev.off(); dev.off()


  
}, future.chunk.size=1 )
})

