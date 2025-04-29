# Chargement des bibliothèques nécessaires
library(dplyr)
library(TraMineR)
library(dplyr)

# Définition des paramètres pour la génération des données
set.seed(123)  # Pour la reproductibilité
n_individus <- 100  # Nombre d'individus
n_annees <- 46  # Nombre d'années (de 1978 à 2023)

# Génération des NIR (identifiants uniques)
NIR <- sprintf("NIR%03d", 1:n_individus)

# Génération des données
data <- data.frame(
  NIR = rep(NIR, each = n_annees),
  annee = rep(1978:2023, times = n_individus),
  salaire = runif(n_individus * n_annees, 1500, 7000),  # Salaire aléatoire entre 1500 et 7000
  sexe = rep(sample(c("Homme", "Femme"), n_individus, replace = TRUE), each = n_annees),
  grade = rep(sample(c("Grade1", "Grade2", "Grade3"), n_individus, replace = TRUE), each = n_annees),
  categorie = rep(sample(c("CatA", "CatB", "CatC"), n_individus, replace = TRUE), each = n_annees),
  statut = rep(sample(c("Titulaire", "Non titulaire"), n_individus, replace = TRUE), each = n_annees),
  age = rep(sample(18:65, n_individus, replace = TRUE), each = n_annees) + rep(0:(n_annees-1), times = n_individus),
  departement = rep(sample(1:95, n_individus, replace = TRUE), each = n_annees),
  region = rep(sample(1:18, n_individus, replace = TRUE), each = n_annees)
)

# Affichage des premières lignes du jeu de données
head(data)

# Supposons que votre base de données soit un dataframe appelé 'data'
# Avec des colonnes 'NIR', 'annee', 'salaire', etc.

# Sélection des individus présents au moins 10 ans
individus_selectionnes <- data %>%
  group_by(NIR) %>%
  filter(n() >= 10) %>%
  ungroup()

# Définition des tranches de salaire
data_with_states <- individus_selectionnes %>%
  mutate(etat = case_when(
    salaire < 2000 ~ "Tranche1",
    salaire >= 2000 & salaire < 4000 ~ "Tranche2",
    salaire >= 4000 & salaire < 6000 ~ "Tranche3",
    salaire >= 6000 ~ "Tranche4"
  ))

# Attribution des états aux individus pour chaque période de temps
seq_wide <- data_with_states %>%
  select(NIR, annee, etat) %>%
  pivot_wider(names_from = annee, values_from = etat) %>%
  arrange(NIR)

# Conversion des trajectoires en séquences
seq_data <- seqdef(seq_wide[, -1], id = seq_wide$NIR, xtstep = 1)

# Calcul des distances entre les séquences
submat <- seqsubm(seq_data, method = "TRATE")

# Calcul de la matrice de distance avec OM
dist_matrix <- seqdist(seq_data, method = "OM", sm = submat, indel = 1)

# Clustering des séquences
clusters <- hclust(as.dist(dist_matrix), method = "ward.D2")

# Visualisation des clusters
plot(clusters) # Dendogramme illisible

# Par exemple, pour 4 groupes
cluster_membership <- cutree(clusters, k = 4)

plot(clusters, labels = FALSE, hang = -1, main = "Dendrogramme des trajectoires")
rect.hclust(clusters, k = 4, border = "red")  # Pour afficher les groupes

# Ajout des groupes aux données
seq_data$cluster <- factor(cluster_membership)
seqIplot(seq_data, group = seq_data$cluster, sortv = "from.start", border = NA)
seqmsplot(seq_data, group = seq_data$cluster, sortv = "from.start")
seqdplot(seq_data, group = seq_data$cluster)

seqdplot(seq_data, group = seq_data$cluster)     # Répartition des états dans le temps
seqmsplot(seq_data, group = seq_data$cluster)    # Trajectoire moyenne
seqIplot(seq_data, group = seq_data$cluster)     # Trajectoires individuelles

seqmsplot(seq_data, group = seq_data$cluster, main = "Trajectoire moyenne par groupe")

## Visualisation par graphiques 
# Liste des centroïdes par groupe
group_medoids <- list()
medoid_indices <- c()

for (k in levels(seq_data$cluster)) {
  # Indices du groupe k
  indices <- which(seq_data$cluster == k)
  
  # Sous-matrice de distances
  sub_dist <- dist_matrix[indices, indices]
  
  # Index de la séquence la plus centrale (somme de distances minimale)
  medoid_index <- indices[which.min(rowSums(as.matrix(sub_dist)))]
  
  # Sauvegarder l’indice
  medoid_indices <- c(medoid_indices, medoid_index)
  group_medoids[[k]] <- seq_data[medoid_index, ]
}

# Afficher les trajectoires centrales
par(mfrow = c(2, 2))  # pour 4 groupes
for (i in seq_along(medoid_indices)) {
  seqplot(seq_data, type = "I", idxs = medoid_indices[i], main = paste("Medoid du groupe", i))
}
par(mfrow = c(1, 1))


## Visualisation en lettres
# Extraire les séquences typiques comme vecteurs
typical_sequences <- lapply(medoid_indices, function(i) {
  as.character(seq_data[i, 1:ncol(seq_data) - 1])  # -1 pour ne pas inclure la colonne "cluster"
})

# Les convertir en chaînes lisibles
typical_labels <- lapply(typical_sequences, function(seq) {
  paste(seq, collapse = " → ")
})

# Afficher les résultats
for (i in seq_along(typical_labels)) {
  cat(paste0("Groupe ", i, " : ", typical_labels[[i]], "\n"))
}

