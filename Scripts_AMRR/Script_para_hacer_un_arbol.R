

### Árbol filogenético

### Descargue y abri como stringset los 10 organismos mas cercanos del blastn

readDNAStringSet("Secuencias_T2/Organismos_cercanos_dino_blast.fasta") -> Organismos_cercanos_dino

### Corroboré que no hubieran caracteres especiales 

alphabetFrequency(Organismos_cercanos_dino)

### Utilicé msa para hacer el alineamiento de las secuencias de DNA

library(msa)

### Alineé las secuencias con clustalW

cercanos_dino_clustalw <- msa(Organismos_cercanos_dino, method = "ClustalW")
cercanos_dino_clustalw

### Convertí el alineamiento a StringSet para poder hacer después la matriz de distancias

dino_alineado <- DNAStringSet(cercanos_dino_clustalw)
dino_alineado

library(DECIPHER)

### Objeto utilizado para la matriz de distancias 

cercanos_dino_clustalw

### Con la función DistanceMatrix hice la matriz de distancias por el método overlap

matriz_distancias_dino <- DistanceMatrix(dino_alineado, method = "overlap")  

View(matriz_distancias_dino)

### Guardé la matriz de distancias en un archivo CSV

write.csv(matriz_distancias_dino,"Resultados_T2/matriz_distancias_dino_T2_P6_Martinez_Melissa")


### La librería ape, disponible en Bioconductor, permite generar árboles filogenéticos

library(ape)

### El árbol lo guarde en formato PDF, en resultados.

pdf("Resultados_T2/arbol_filogenetico_Martinez_Melissa_T2_P6.pdf", width = 8, height = 5)

### Creé un árbol filogenético utilizando el método Neighbor-Joining
### a partir de la matriz de distancias

arbol_filogenetico_dino <- nj(matriz_distancias_dino)

### Esta función nos ayuda a personalizar el árbol.

plot( arbol_filogenetico_dino, 
      cex = 0.6, # tamaño de letra
      type = "cladogram", # tipo de arbol
      tip.color = "darkblue", # color de letra
      edge.color = "pink", # color de ramas    
      edge.width = 3) # tamaño de ramas


dev.off()