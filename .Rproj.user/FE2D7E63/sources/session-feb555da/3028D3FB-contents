library(Biostrings)

### Árbol filogenético

### Abrí mi documento fasta 
readDNAStringSet("Secuencias_AMMR/todos_los_virus.fasta") -> virus

### Corroboré que no hubieran caracteres especiales 

alphabetFrequency(virus)

### Utilicé msa para hacer el alineamiento de las secuencias de DNA

library(msa)

### Alineé las secuencias con clustalW

virus_clustalw <- msa(virus, method = "ClustalW")
virus_clustalw

### Alineé las secuencias con muscle

virus_muscle <- msa(virus, method = "Muscle")
virus_muscle

### Convertí el alineamiento en un formato que pudiera utilizar para hacer una matriz de distancias

virus_alineado_clustalw <- msaConvert(virus_clustalw, type = "seqinr::alignment")

virus_alineado_muscle <- msaConvert(virus_muscle, type = "seqinr::alignment")

### Libreria para hacer matriz de distancias

library(seqinr)

# Calculo la matriz de distancias 
distancias_virus_clustalw <- dist.alignment(virus_alineado_clustalw, "identity")

distancias_virus_muscle <- dist.alignment(virus_alineado_muscle, "identity")

matriz_dist_clustalw <- as.matrix(distancias_virus_clustalw)
matriz_dist_clustalw

matriz_dist_muscle <- as.matrix(distancias_virus_muscle)
matriz_dist_muscle

### La librería ape, disponible en Bioconductor, permite generar árboles filogenéticos

library(ape)

### El árbol lo guarde en formato PDF, en resultados.

pdf("Resultados_AMMR/arbol_filogenetico_virus_clustalw.pdf", width = 8, height = 5)

### Creé un árbol filogenético utilizando el método Neighbor-Joining, primero con
### el alineamiento de clustalw a partir de la matriz de distancias

arbol_filogenetico_virus <- nj(matriz_dist_clustalw)

### Esta función nos ayuda a personalizar el árbol.

plot( arbol_filogenetico_virus, 
      cex = 0.6, # tamaño de letra
      type = "cladogram", # tipo de arbol
      tip.color = "darkblue", # color de letra
      edge.color = "pink", # color de ramas    
      edge.width = 3) # tamaño de ramas


dev.off()

pdf("Resultados_AMMR/arbol_filogenetico_virus_muscle.pdf", width = 8, height = 5)

### Creé un árbol filogenético utilizando el método Neighbor-Joining ahoara 
### para el alineamiento de muscle, a partir de la matriz de distancias

arbol_filogenetico_virus2 <- nj(matriz_dist_muscle)

### Esta función nos ayuda a personalizar el árbol.

plot( arbol_filogenetico_virus2, 
      cex = 0.6, # tamaño de letra
      type = "cladogram", # tipo de arbol
      tip.color = "darkblue", # color de letra
      edge.color = "pink", # color de ramas    
      edge.width = 3) # tamaño de ramas


dev.off()

## Nuevo arbol usando ggtree

arbol_filogenetico_virus_clustalw <- nj(matriz_dist_clustalw)

library(ggplot2)  
library(ggtree)

## Arbol con ggtree "normal"

ggtree_plot <- ggtree(arbol_filogenetico_virus_clustalw) +
  geom_tiplab(size = 5, color = "darkblue") +  
  geom_tree(color = "green", size = 3) +   
  xlim(0,0.5)

## Guardar en PDF
ggsave("Resultados_AMMR/arbol_filogenetico_virus_clustalw_ggtree.pdf",
       plot = ggtree_plot, width = 40, height = 15)

## Arbol circular

ggtree_plot <- ggtree(arbol_filogenetico_virus_clustalw, layout = "circular") +
  geom_tiplab(size = 3, color = "purple") +
  geom_tree(color = "pink", size = 1) +
  xlim(0,0.5)

## Guardar el arbol circular

ggsave("Resultados_AMMR/arbol_filogenetico_virus_clustalw_ggtree_circular.pdf", 
       plot = ggtree_plot, width = 20, height = 20) 

