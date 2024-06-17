## Cargar librerias necesarias para el análisis ####
# Lista de librerias necesarias
packages_list<-list("magrittr", "dplyr", "plyr", "this.path", "ggplot2", "raster", "terra", "sf", "ggspatial",
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
print(dir_work)

## Cargar area de estudio
studyArea<- terra::rast(file.path(dir_work, "StudyArea", "studyAreaCAR.tif")) %>% setNames("StudyArea")

## Cargar opciones area M ####
areaM_preliminar<-  terra::rast(file.path(dir_work, "AreaM", "areaM_buffer1grado.tif")) %>% setNames("areaM_buffer1grado")
areaM<-  terra::rast(file.path(dir_work, "AreaM", "areaM_convexo_OrobiomaAndinoAltoandinocordilleraoriental.tif")) %>% setNames("AreaM")

## Cargar puntos registros
data_agras<- read.csv(file.path(dir_work, "dataOccurrences" ,  "BDmodelo_CAR.csv"))
data_agras_points <- data_agras %>% st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs= 4326) 

## Establecer covariables ####
dir_layers<- file.path(dir_work, "variables_AreaM") # folder donde se almacenan todas las covariables
explore_layers <- list.files(dir_layers, recursive = F)
print(explore_layers)

## Cargar variables ####

env_layers<- list.files( dir_layers, "\\.asc$|\\.tif$", recursive = T, full.names = T)
envMstack<- terra::rast(env_layers)
names(envMstack)<-  gsub(",", ".", names(envMstack)); names(envMstack)<-  gsub("-", "",  names(envMstack)) # ajustar nombres

# reescribir covariables a dormato asc para optimizar el analisis
folder_vars_areaM<- file.path(dir_work, "Variables_AreaM_adjust/Set_1") # folder donde se reescriben

# generar folder de variables ajustadas
dir.create(folder_vars_areaM, showWarnings = F, recursive = T)
files_envM<- list.files( folder_vars_areaM, "\\.asc$", recursive = T, full.names = F) %>% tools::file_path_sans_ext()

# funcion de ajuste
ndec<- function(x) {dec <- nchar(strsplit(as.character(x), "\\.")[[1]][2]);
if (is.na(dec)) dec <- 0
return(as.numeric(dec))}

# loop ajuste por variable
for (i in 1:nlyr(envMstack)) {  name_layer<- names(envMstack)[i] ;
if(!(name_layer %in% files_envM)){ # si la variable se ajusto previamente se salta este paso
  
  envMi <- envMstack[[i]]
  
  value <- values(envMi, na.rm = TRUE) %>% max()
  y <- ndec(x = value)
  
  if (y == 0) { datTyp <- "INT2S"; decinum <- 0 }
  if (y >= 1) { datTyp <- "FLT4S"; decinum <- 3 }
  
  writeRaster(x = if(is.factor(envMi)){envMi}else{round(envMi, digits = decinum)},
              filename = file.path(folder_vars_areaM, paste0( names(envMstack[[i]]), ".asc")), overwrite = T, NAflag = -9999, datatype = datTyp)
}   }

# Cargar stac de datos ajustados
env.Mfiles <- list.files(folder_vars_areaM, ".asc$", recursive = T, full.names = T)
env.M <- terra::rast(env.Mfiles)
print(env.M)

# umbral de correlacion para pruebas de multicolinealidad 
cor_threshold<- 0.65

# obtener covariables por ocurrencia ####

## Crear grilla base de analisis
raster_base<- terra::rast(areaM_preliminar)
raster_ids<- raster_base %>%  terra::setValues(seq(ncell(.)))
occurrences_ids<- terra::extract(raster_ids, data_agras_points, ID= F) %>% {.[!duplicated(.[,1]),]}


## data ocurrencias variables
data_occs<- env.M[occurrences_ids]  %>% 
  select_if(~ !all(is.na(.))) %>%  # eliminar variables no informativas (todo NA)
  mutate_if(~ all(. %in% c(0, 1), na.rm = TRUE), ~ as.factor(as.character(.))) %>% # convertir variables binarias en factor
  select_if(~ !is.factor(.) || (is.factor(.) && nlevels(.) > 1)) %>% # remover factores no informativos (1 solo nivel)
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>% # Inferir valores faltantes - media numerica
  mutate(across(where(is.factor), ~ {mode_val <- names(sort(table(.), decreasing = TRUE))[1];replace_na(., mode_val)})) %>% # Inferir valores faltantes - nivel mas comun del factor
  {cbind(setNames(as.data.frame(terra::xyFromCell(raster_ids, occurrences_ids)), c("longitude", "latitude")), .)} # Asignar coordenadas

# Ajustar stack variables informativas
envMstack<- env.M[[ names(env.M) %>% {.[.%in% names(data_occs)]}  ]] # filtrar por stack de variables informativas





# numero de datos aleatorios - pseudo ausencias 
Max.Bg<- 10000

# variables que deben mantenerse en pruebas de multicolinealidad ####
vars_predilection<-  c("topographic_earthenv_Elevation", "worldclim_bio12_Annual_Precipitation", "CLC_CLC2_Arbustales_rep_apot0.5km") %>% 
  {data.frame(Var= ., sort_pred=  seq_along(.))}

# Definición parametros prueba modelos 
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

# tune list
tune_list<- list(fc = unique(models_input$fc), rm = as.numeric(unique(models_input$rm)))

#  proc
proc_userval_prel<- function(vars) {
  proc <- kuenm::kuenm_proc(vars$occs.val.pred, c(vars$bg.train.pred, vars$bg.val.pred))
  out <- data.frame(
    proc_auc_ratio = proc$pROC_summary[1],
    proc_pval = proc$pROC_summary[2], row.names = NULL
  )
  return(out)
}



## Cargar mapa de expertos ####
mapaExpertos<- terra::rast(file.path(dir_work, "ExpertsMap", "mapaExpertos.tif")) %>% setNames("mapaExperto")
mapaExpertos_pol<- mapaExpertos %>% terra::as.polygons() %>% sf::st_as_sf() %>% dplyr::mutate(level_map= c("Posible\nAusencia", "Posible\nPresencia"))

plot_mapaExpertos<- ggplot()+
  annotation_map_tile(type = "cartolight", zoom = 7) +
  geom_sf(data = st_as_sf(as.polygons(areaM)), aes(color = "orange"), fill = NA, size= 0.1 )+
  scale_color_identity(name = "", guide = "legend",  labels = c( "Área M"),breaks = c("orange"))+
  ggnewscale::new_scale_fill()+
  geom_sf(data = mapaExpertos_pol, aes(fill = level_map), alpha = 0.5, color = NA)+
  scale_fill_manual("Mapa de\nexpertos", values = setNames(c("purple", "darkgreen"), c("Posible\nAusencia", "Posible\nPresencia")))+
  theme(text = element_text(size = 8))




# Calibracion del modelo - Seleccion iterativa de variables  ####

# definir bg - background data - pseduo ausencias 
ids_background_prel<- cells(areaM_preliminar) %>% {.[!. %in% occurrences_ids ]}
sample_background_prel<- sample(x = ids_background_prel, size = Max.Bg, replace = F)
data_Sbg_prel <- envMstack %>% {.[sample_background_prel]} %>% # cargar celdas
  dplyr::select(intersect(names(.), names(data_occs))) %>% mutate(across(everything(), ~ { # ajustar str variables acorde ocurrencias
    if (is.factor(data_occs[[cur_column()]])) { factor(., levels = levels(data_occs[[cur_column()]]))
      } else if (is.numeric(data_occs[[cur_column()]])) { as.numeric(.)
        } else {.} })) %>% 
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>% # Inferir valores faltantes - media numerica
  mutate(across(where(is.factor), ~ {mode_val <- names(sort(table(.), decreasing = TRUE))[1];replace_na(., mode_val)})) %>% # Inferir valores faltantes - nivel mas comun del factor
  {cbind(setNames(as.data.frame(terra::xyFromCell(raster_ids, sample_background_prel)), c("longitude", "latitude")), .)} # Asignar coordenadas

# organizar data test multicol. Requiere binarias como numericas
vars_data_cor_prel<- envMstack[terra::cells(areaM_preliminar)]  %>%
  dplyr::select(intersect(names(.), names(data_occs))) %>% mutate(across(everything(), ~ { # ajustar str variables acorde ocurrencias
    if (is.factor(data_occs[[cur_column()]])) { factor(., levels = levels(data_occs[[cur_column()]]))
    } else if (is.numeric(data_occs[[cur_column()]])) { as.numeric(.)
    } else {.} })) %>% 
  mutate_if(is.factor, ~ as.numeric(as.character(.))) %>% # convertir factores en numericas
  mutate_if(~ !is.numeric(.) && !is.factor(.), as.numeric)
  
vars_test_prel <- names(vars_data_cor_prel) # objeto dinamico de variables a probar/ se actualizara dentro del loop de prueba

## loop iterativo multicolinealidad - importancia de variables
list_hclust_model_prel<-  list()     # registro iterativo de cada ejecucion
vars_test_prel_vif_model<- vars_test_prel

for(j in seq_along(vars_test_prel_vif_model)){ print(paste0("j_", j))
  
  ## Test de multicolinealidad ####
  for(v in seq_along(vars_test_prel_vif_model)){ print(paste0("v_", v))
    
    vars_data_test_model_prel<- dplyr::select(vars_data_cor_prel, vars_test_prel_vif_model)
    
    ## Prueba de correlacion ####
    cordata_model_prelR_prel<- cor(vars_data_test_model_prel, method = "spearman") %>% as.data.frame.matrix()  %>% set_rownames(names(.))
    NACol<- names(which(rowSums(is.na(cordata_model_prelR_prel)) > (ncol(cordata_model_prelR_prel)/2) ))
    cordata_model_prel<- cordata_model_prelR_prel %>% {.[!names(.) %in% NACol,]} %>% {.[,!colnames(.) %in% NACol]}; 
    cordata_model_prel[is.na(cordata_model_prel)]<-0
    
    corhclust_model_prel <- hclust(as.dist(1-abs(cordata_model_prel)))
    group_covars_model_prel<- cutree(corhclust_model_prel, h = 1-cor_threshold) %>% as.data.frame %>% rownames_to_column("Var") %>%
      data.table::setnames(".", "group")
    
    ## Estimacion VIF ####
    vif_data_model_prel <- ginv(as.matrix(cordata_model_prel)) %>% diag() %>% setNames(names(cordata_model_prel)) %>% as.data.frame() %>% rownames_to_column("Var") %>% 
      setNames(c("Var", "VIF")) %>% list(group_covars_model_prel,vars_predilection) %>% plyr::join_all() %>% arrange(group, sort_pred, VIF)
    
    list_hclust_model_prel[[length(list_hclust_model_prel)+1]]<- list(corhclust=corhclust_model_prel, vif_data= vif_data_model_prel, cor_data= cordata_model_prel)
    
    check_duplicated_prel<- any(duplicated(group_covars_model_prel$group))
    
    if( check_duplicated_prel ){
      
      vars_test_prel_vif_model<- vif_data_model_prel %>% dplyr::filter(!duplicated(group))  %>% {.$Var}
      
    } else {
      break
    }
  }
  
  
  ## Test de calibracion del modelo  ####
  data_occs_test_prel<- data_occs %>% dplyr::select(c("longitude", "latitude", vars_test_prel_vif_model))  
  Sbg_test_prel<- data_Sbg_prel %>% dplyr::select(c("longitude", "latitude", vars_test_prel_vif_model)) 

  
  eval_block_prel <- ENMevaluate(
    occs = data_occs_test_prel , bg = Sbg_test_prel, partitions = c("block"),
    tune.args = tune_list, algorithm = "maxent.jar",
    doClamp = T, user.eval = proc_userval_prel, parallel = T, numCores= 25
  )
  
  future::plan(sequential); gc()
  length_models<- length(eval_block_prel@variable.importance)
  
  ### Importancia de variables entre modelos ####
  importance_data_prel <- eval_block_prel@variable.importance %>% plyr::rbind.fill() %>% dplyr::group_by(variable) %>% 
    dplyr::summarise(
      percent.contribution = sum(percent.contribution)/length_models,
      permutation.importance = sum(permutation.importance)/length_models,
    ) %>% arrange(-permutation.importance) %>% dplyr::rename(Var=variable)
  
  list_hclust_model_prel[[length(list_hclust_model_prel)]]$importance<- importance_data_prel
  list_hclust_model_prel[[length(list_hclust_model_prel)]]$vif_data<- list(importance_data_prel, vif_data_model_prel) %>% plyr::join_all()
  
  if( !any(importance_data_prel$permutation.importance<1)  ){
    break
  } else {
    ### Eliminacion de variables no informativas 1% ####
    vars_test_prel_vif_model<- importance_data_prel %>% as.data.frame() %>%  arrange(permutation.importance) %>% dplyr::filter(permutation.importance>=1) %>% {.$Var}
  }
}

  list_selectvars_model_prel <- pblapply(list_hclust_model_prel, function(y) {
    
    # Convertir agrupamiento gerarquico en dendograma
    cordend_data <- dendro_data(as.dendrogram(y$corhclust))
    
    # Organizar datos vif por variables
    vif_data_h <- y$vif_data
    
    # Verificar si la prueba estimo modelo de importancia por validacion cruzada o fue una prueba preliminar de multicolinealidad
    check_permut <- "permutation.importance" %in% names(vif_data_h)
    
    # Crear una tabla de variables con información del dendrograma.
    var_table <- with(cordend_data$labels, data.frame(y_center = x, y_min = x - 0.5, y_max = x + 0.5, Var = as.character(label), height = 1)) %>% 
      list(vif_data_h) %>% plyr::join_all() %>%  arrange(group, VIF) %>% 
      mutate(Var = if_else(!duplicated(group), paste("*", Var , sep = ""), Var),
             colVar = if_else(!duplicated(group), "red", "black")) %>%  # resaltar colores variables que superan multicol
      arrange(y_center) %>% 
      mutate(change_flag = if_else(group != lag(group, default = first(group)), 1, 0),
             intercalated = cumsum(change_flag) %% 2) %>% 
      dplyr::mutate(col = if_else(intercalated == 1, "#EBEBEB", "white"))
    
    # Si se estimo importancia, actualizar la tabla de variables.
    if (check_permut) {
      var_table <- var_table %>% 
        mutate(Var = if_else(permutation.importance >= 1, paste("*", Var, sep = ""), Var),
               colVar = if_else(permutation.importance >= 1, "red", "black"))  # resaltar colores variables que superan importancia
    }
    
    # Crear datos de segmentos para el dendrograma.
    segment_data <- with(segment(cordend_data), data.frame(x = y, y = x, xend = yend, yend = xend, cor = 1 - yend))
    
    # Crear el plot del dendrograma
    ggdendroPlot <- ggplot() +
      annotate("rect", xmin = -0.05, xmax = 1.04, fill = var_table$col, ymin = var_table$y_min, ymax = var_table$y_max, alpha = 0.9) +
      geom_segment(data = segment_data, aes(x = 1 - x, y = y, xend = 1 - xend, yend = yend, label = cor), size = 0.3) +
      scale_y_continuous(breaks = var_table$y_center, labels = var_table$Var,
                         sec.axis = sec_axis(~ ., name = "VIF", breaks = var_table$y_center, labels = round(var_table$VIF, 1))) +
      coord_cartesian(expand = FALSE) +
      labs(x = "Correlation", y = "Variables") +
      geom_vline(xintercept = cor_threshold, linetype = "dashed", col = "red") +
      theme(legend.position = "bottom", legend.key.width = unit(50, 'pt'),
            plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
            panel.grid.major = element_line(color = "gray"),
            axis.ticks.length = unit(0.3, "mm"),
            text = element_text(size = 10),
            axis.text.y = element_text(colour = var_table$colVar))
    
    # Si se estimo importancia, agregar etiquetas de importancia y permutacion
    if (check_permut) {
      labels_importance <- var_table %>% 
        dplyr::mutate(label_imp = paste0("[", round(percent.contribution), ";", round(permutation.importance), "]"), 
                      label_imp_col = if_else(permutation.importance >= 1, "blue", "black"))
      
      ggdendroPlot <- ggdendroPlot +
        geom_label(data = labels_importance, aes(y = y_center, x = Inf, label = label_imp), 
                   color = labels_importance$label_imp_col, size = 3, hjust = 1.1, fill = var_table$col, label.size = 0)
    }
    
    # Agregar el gráfico en la lista.
    y$ggdendroPlot <- ggdendroPlot
    
    # Devolver el elemento modificado de la lista.
    y
  })



list_selectvars_model_prel[[length(list_hclust_model_prel)]]$corhclust
ggsave(file.path(dirname(dir_work), "README_figures", paste0("compiled_dendrogram_selectvars_model_prel", ".png")) ,
       list_selectvars_model_prel[[length(list_selectvars_model_prel)]][["ggdendroPlot"]],
      width= 10, height = 4
       )




saveRDS(eval_block_prel, "eval_block_prel.rds")

saveRDS(compiled_dendrogram_selectvars_model_prel, "compiled_dendrogram_selectvars_model_prel.rds")


# Obtener compilacion de dendogramas de todas las listas
compiled_dendrogram_selectvars_model_prel <- purrr::map(list_selectvars_model_prel, "ggdendroPlot")

# Exportar pdf compilado
pdf(file.path(dirname(dir_work), "README_figures", paste0("compiled_dendrogram_selectvars_model_prel", ".pdf")))
for(plot in compiled_dendrogram_selectvars_model_prel) {
  print(plot)
}

tryCatch({dev.off(); dev.off()}, error= function(e) NULL)



## Figura Importancia de variables ####

# Extraer info del modelo
eval_importance_model<- eval_block_prel@variable.importance

# Organizar información y estimar la métrica de importancia entre modelos
length_models<- length(eval_block_prel@variable.importance)
summ_importance_data_prel <- eval_block_prel@variable.importance %>% plyr::rbind.fill() %>% dplyr::group_by(variable) %>% 
  dplyr::summarise(
    percent.contribution = sum(percent.contribution)/length_models,
    permutation.importance = sum(permutation.importance)/length_models,
  ) %>%   dplyr::rowwise() %>% dplyr::mutate(ranking= mean(c(percent.contribution, permutation.importance))) %>% 
  arrange(ranking) %>% 
  as.data.frame() %>% dplyr::mutate(variable= factor(variable, levels = unique(.$variable)))

# Ajustar los datos para el gráfico
dataplot_impvars_model_prel<-  summ_importance_data_prel %>% 
  pivot_longer(cols = c(percent.contribution, permutation.importance), 
               names_to = "tipo", values_to = "aporte") %>% 
  dplyr::mutate(tipo=factor(tipo, levels= rev(c("percent.contribution", "permutation.importance"))))

# Gráfico de Importancia de Variables
plot_impvars_model_prel <-     ggplot(dataplot_impvars_model_prel, aes(x = variable)) +
  geom_bar(aes(y = aporte, fill = tipo), stat = "identity", alpha= 0.5)+
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100))+
  coord_flip() + 
  labs(x = "Variables", y = "Porcentaje", fill = "") +
  scale_fill_manual(values = setNames(c("#0000FE", "#00ADAC"), 
                                      c("percent.contribution", "permutation.importance")),
                    labels = setNames(c("Contribución", "Pérdida\npermutación"), 
                                      c("percent.contribution", "permutation.importance"))
  ) + theme_minimal()


  
ggsave(file.path(dirname(dir_work), "README_figures", paste0("plot_impvars_model_prel", ".png")) ,
       plot_impvars_model_prel,
       width= 10, height = 4
)





# En Bioma andino M ##############################################################################

# Calibracion del modelo - Seleccion iterativa de variables  ####

# definir bg - background data - pseduo ausencias 
ids_background<- cells(areaM) %>% {.[!. %in% occurrences_ids ]}
sample_background<- sample(x = ids_background, size = Max.Bg, replace = F)
data_Sbg <- envMstack %>% {.[sample_background]} %>% # cargar celdas
  dplyr::select(intersect(names(.), names(data_occs))) %>% mutate(across(everything(), ~ { # ajustar str variables acorde ocurrencias
    if (is.factor(data_occs[[cur_column()]])) { factor(., levels = levels(data_occs[[cur_column()]]))
    } else if (is.numeric(data_occs[[cur_column()]])) { as.numeric(.)
    } else {.} })) %>% 
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>% # Inferir valores faltantes - media numerica
  mutate(across(where(is.factor), ~ {mode_val <- names(sort(table(.), decreasing = TRUE))[1];replace_na(., mode_val)})) %>% # Inferir valores faltantes - nivel mas comun del factor
  {cbind(setNames(as.data.frame(terra::xyFromCell(raster_ids, sample_background)), c("longitude", "latitude")), .)} # Asignar coordenadas

# organizar data test multicol. Requiere binarias como numericas
vars_data_cor<- envMstack[terra::cells(areaM)]  %>%
  dplyr::select(intersect(names(.), names(data_occs))) %>% mutate(across(everything(), ~ { # ajustar str variables acorde ocurrencias
    if (is.factor(data_occs[[cur_column()]])) { factor(., levels = levels(data_occs[[cur_column()]]))
    } else if (is.numeric(data_occs[[cur_column()]])) { as.numeric(.)
    } else {.} })) %>% 
  mutate_if(is.factor, ~ as.numeric(as.character(.))) %>% # convertir factores en numericas
  mutate_if(~ !is.numeric(.) && !is.factor(.), as.numeric)

vars_test <- names(vars_data_cor) # objeto dinamico de variables a probar/ se actualizara dentro del loop de prueba

## loop iterativo multicolinealidad - importancia de variables
list_hclust_model<-  list()     # registro iterativo de cada ejecucion
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
  data_occs_test<- data_occs %>% dplyr::select(c("longitude", "latitude", vars_test_vif_model))  
  Sbg_test<- data_Sbg %>% dplyr::select(c("longitude", "latitude", vars_test_vif_model)) 
  
  
  eval_block <- ENMevaluate(
    occs = data_occs_test , bg = Sbg_test, partitions = c("block"),
    tune.args = tune_list, algorithm = "maxent.jar",
    doClamp = T, user.eval = proc_userval, parallel = T, numCores= 25
  )
  
  future::plan(sequential); gc()
  length_models<- length(eval_block@variable.importance)
  
  ### Importancia de variables entre modelos ####
  importance_data <- eval_block@variable.importance %>% plyr::rbind.fill() %>% dplyr::group_by(variable) %>% 
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

list_selectvars_model <- pblapply(list_hclust_model, function(y) {
  
  # Convertir agrupamiento gerarquico en dendograma
  cordend_data <- dendro_data(as.dendrogram(y$corhclust))
  
  # Organizar datos vif por variables
  vif_data_h <- y$vif_data
  
  # Verificar si la prueba estimo modelo de importancia por validacion cruzada o fue una prueba preliminar de multicolinealidad
  check_permut <- "permutation.importance" %in% names(vif_data_h)
  
  # Crear una tabla de variables con información del dendrograma.
  var_table <- with(cordend_data$labels, data.frame(y_center = x, y_min = x - 0.5, y_max = x + 0.5, Var = as.character(label), height = 1)) %>% 
    list(vif_data_h) %>% plyr::join_all() %>%  arrange(group, VIF) %>% 
    mutate(Var = if_else(!duplicated(group), paste("*", Var , sep = ""), Var),
           colVar = if_else(!duplicated(group), "red", "black")) %>%  # resaltar colores variables que superan multicol
    arrange(y_center) %>% 
    mutate(change_flag = if_else(group != lag(group, default = first(group)), 1, 0),
           intercalated = cumsum(change_flag) %% 2) %>% 
    dplyr::mutate(col = if_else(intercalated == 1, "#EBEBEB", "white"))
  
  # Si se estimo importancia, actualizar la tabla de variables.
  if (check_permut) {
    var_table <- var_table %>% 
      mutate(Var = if_else(permutation.importance >= 1, paste("*", Var, sep = ""), Var),
             colVar = if_else(permutation.importance >= 1, "red", "black"))  # resaltar colores variables que superan importancia
  }
  
  # Crear datos de segmentos para el dendrograma.
  segment_data <- with(segment(cordend_data), data.frame(x = y, y = x, xend = yend, yend = xend, cor = 1 - yend))
  
  # Crear el plot del dendrograma
  ggdendroPlot <- ggplot() +
    annotate("rect", xmin = -0.05, xmax = 1.04, fill = var_table$col, ymin = var_table$y_min, ymax = var_table$y_max, alpha = 0.9) +
    geom_segment(data = segment_data, aes(x = 1 - x, y = y, xend = 1 - xend, yend = yend, label = cor), size = 0.3) +
    scale_y_continuous(breaks = var_table$y_center, labels = var_table$Var,
                       sec.axis = sec_axis(~ ., name = "VIF", breaks = var_table$y_center, labels = round(var_table$VIF, 1))) +
    coord_cartesian(expand = FALSE) +
    labs(x = "Correlation", y = "Variables") +
    geom_vline(xintercept = cor_threshold, linetype = "dashed", col = "red") +
    theme(legend.position = "bottom", legend.key.width = unit(50, 'pt'),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
          panel.grid.major = element_line(color = "gray"),
          axis.ticks.length = unit(0.3, "mm"),
          text = element_text(size = 10),
          axis.text.y = element_text(colour = var_table$colVar))
  
  # Si se estimo importancia, agregar etiquetas de importancia y permutacion
  if (check_permut) {
    labels_importance <- var_table %>% 
      dplyr::mutate(label_imp = paste0("[", round(percent.contribution), ";", round(permutation.importance), "]"), 
                    label_imp_col = if_else(permutation.importance >= 1, "blue", "black"))
    
    ggdendroPlot <- ggdendroPlot +
      geom_label(data = labels_importance, aes(y = y_center, x = Inf, label = label_imp), 
                 color = labels_importance$label_imp_col, size = 3, hjust = 1.1, fill = var_table$col, label.size = 0)
  }
  
  # Agregar el gráfico en la lista.
  y$ggdendroPlot <- ggdendroPlot
  
  # Devolver el elemento modificado de la lista.
  y
})



ggsave(file.path(dirname(dir_work), "README_figures", paste0("compiled_dendrogram_selectvars_model", ".png")) ,
       list_selectvars_model[[length(list_selectvars_model)]][["ggdendroPlot"]],
       width= 10, height = 4
)




# Obtener compilacion de dendogramas de todas las listas
compiled_dendrogram_selectvars_model <- purrr::map(list_selectvars_model, "ggdendroPlot")

# Exportar pdf compilado
pdf(file.path(dirname(dir_work), "README_figures", paste0("compiled_dendrogram_selectvars_model", ".pdf")))
for(plot in compiled_dendrogram_selectvars_model) {
  print(plot)
}

tryCatch({dev.off(); dev.off()}, error= function(e) NULL)


eval_block<- eval_jack
## Figura Importancia de variables ####

# Extraer info del modelo
eval_importance_model<- eval_block@variable.importance

# Organizar información y estimar la métrica de importancia entre modelos
length_models<- length(eval_block@variable.importance)
summ_importance_data <- eval_block@variable.importance %>% plyr::rbind.fill() %>% dplyr::group_by(variable) %>% 
  dplyr::summarise(
    percent.contribution = sum(percent.contribution)/length_models,
    permutation.importance = sum(permutation.importance)/length_models,
  ) %>%   dplyr::rowwise() %>% dplyr::mutate(ranking= mean(c(percent.contribution, permutation.importance))) %>% 
  arrange(ranking) %>% 
  as.data.frame() %>% dplyr::mutate(variable= factor(variable, levels = unique(.$variable)))

# Ajustar los datos para el gráfico
dataplot_impvars_model<-  summ_importance_data %>% 
  pivot_longer(cols = c(percent.contribution, permutation.importance), 
               names_to = "tipo", values_to = "aporte") %>% 
  dplyr::mutate(tipo=factor(tipo, levels= rev(c("percent.contribution", "permutation.importance"))))

# Gráfico de Importancia de Variables
plot_impvars_model <-     ggplot(dataplot_impvars_model, aes(x = variable)) +
  geom_bar(aes(y = aporte, fill = tipo), stat = "identity", alpha= 0.5)+
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100))+
  coord_flip() + 
  labs(x = "Variables", y = "Porcentaje", fill = "") +
  scale_fill_manual(values = setNames(c("#0000FE", "#00ADAC"), 
                                      c("percent.contribution", "permutation.importance")),
                    labels = setNames(c("Contribución", "Pérdida\npermutación"), 
                                      c("percent.contribution", "permutation.importance"))
  ) + theme_minimal()



ggsave(file.path(dirname(dir_work), "README_figures", paste0("plot_impvars_model", ".png")) ,
       plot_impvars_model,
       width= 10, height = 4
)













