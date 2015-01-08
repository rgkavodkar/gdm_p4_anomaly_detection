
# read the parameters from the command line
args = commandArgs(TRUE)
# input folder path
input_folder_path = args[1]
# feature type
feature_type = as.integer(args[2])
# the size of the sliding window
window_width = as.integer(args[3])
# value of W'
w_dash = as.integer(args[4])
# output folder path
output_folder_path = args[5]


if(feature_type == "3") {
  library(sna)
  library(igraph)
} else {
  library(igraph)
  library(sna)
  detach("package:sna", unload=TRUE)
}

# Function to return the degrees of the nodes for a given day
get_degrees = function(edges, no_vertices) {
  no_vertices = no_vertices - 1
  day_graph = graph.data.frame(as.matrix(edges), directed=FALSE, vertices=0:no_vertices)
  # gives the degree matrix of the graph
  degree_matrix = degree(day_graph)
  degree_matrix = t(as.matrix(degree_matrix))
  degree_matrix
}

# Function to return the clustering coefficient of the nodes for a given day
get_clustering_coefficient = function(edges, no_vertices) {
  no_vertices = no_vertices - 1
  day_graph = graph.data.frame(as.matrix(edges), directed=FALSE, vertices=0:no_vertices)
  # measure the clustering coefficient
  cc_matrix = transitivity(day_graph, type="local")
  cc_matrix = t(as.matrix(cc_matrix))
  cc_matrix[is.na(cc_matrix)] = 0
  cc_matrix
}

# Function to return the number of nodes in the egonet of the nodes for a given day
get_no_edges_egonet = function(edges, no_vertices) {
  day_graph = graph.data.frame(as.matrix(edges), directed=FALSE, vertices=0:no_vertices)
  adj_matrix = get.adjacency(day_graph)
  adj_matrix = data.matrix(adj_matrix, rownames.force = NA)
  # calculate the egonet of the matrix
  egonet = ego.extract(adj_matrix, neighborhood="combined")
  en_matrix =  matrix(0, ncol=no_vertices,nrow = 1)
  # extract the edge count in the egonet
  for (i in 1:no_vertices){
    index = toString(i - 1)
    egonet_matrix = data.matrix(egonet[[index]])
    en_matrix[1,i] = length(which(egonet_matrix == 1))/2
  }
  en_matrix
  
}

# calculates the moving average required for computing the threshold
calculate_moving_average <- function(values) {
  sum = 0
  len = length(values)
  for (i in 2:len){
    # calculated as abs(a - b) + abs(b - c) + abs(c - d) .. 
    sum = sum + abs( (values[i-1]) - (values[i]) )
  }  
  average = sum/(len -1)
  average
}

# for a given matrix, returns the correlation matrix
get_correlation_matrix <- function(data_matrix) {
  #calculate the pearson correlation matrix
  correlation_matrix = cor(data_matrix, use = "pairwise.complete.obs", method = "pearson")
  correlation_matrix[is.na(correlation_matrix)] = 1
  return(as.matrix(correlation_matrix))
}

# get all the files in the given directory
files = list.files(input_folder_path)

edges = data.frame()
day_data = list()
day_counter = 0

# This snippet of code reads the edges from the time series file and store them into
# a data frame along with a third attr that describes the day number
for(file in files) {
  
  day_counter = day_counter + 1
  
  index = which(strsplit(file, "")[[1]]=="_")
  day = as.integer(substr(file, 1, index[1]-1))
  data=read.table(paste(input_folder_path, file, sep = "/"),header=F, sep=" ")
  no_vertices = data[1,1]
  data = data[-1, ]
  
  # Features and their IDs
  # 1. Degree of the node
  # 2. Clustering Coefficient
  # 3. Number of edges in the node's egonet
  if(feature_type == 1) {
    day_data[[day_counter]] = get_degrees(data, no_vertices)
  } else if(feature_type == 2) {
    day_data[[day_counter]] = get_clustering_coefficient(data, no_vertices)
  } else if(feature_type == 3) {
    day_data[[day_counter]] = get_no_edges_egonet(data, no_vertices)
  }
}

window_data = matrix(0, ncol = no_vertices, nrow = 6)

# initializing the window_data to hold the first 6, so that
# in the coming loop we could just add one to complete the 7 day window
temp = window_width-1
for(i in 1:temp) {
  window_data[i, ] = day_data[[i]]
}

count_correlation_matrix = day_counter - (window_width - 1)
principal_eig_vectors = matrix(0, nrow = no_vertices, ncol = count_correlation_matrix)

current_corr_matrix_index = 1
i = 7
for (i in window_width:day_counter){
  # add the (window_width)th day to the window to complete it
  window_last_day = day_data[[i]]
  window_data = rbind(window_data, window_last_day)
  
  # get the correlation matrix
  correlation_mat = get_correlation_matrix(window_data)
  
  eig_vectors = eigen(correlation_mat)
  eig_vectors = as.matrix(eig_vectors$vectors)

  # consider only the first eigen vector
  principal_eig_vectors[, current_corr_matrix_index] = eig_vectors[, 1]
  current_corr_matrix_index = current_corr_matrix_index + 1

  # remove the first entry from the window to accomodate another in the next loop
  window_data = window_data[-1, ]
}

z = matrix(0, nrow = 1, ncol = count_correlation_matrix)

for (i in (w_dash + 1):count_correlation_matrix){
  # get the previous W' vectors corresponding to r(t-1)
  x = i - 5
  y = i - 1
  
  # ave of last W' principal eigen vectors
  avg_eig_vector = rowMeans(principal_eig_vectors[, x:y])
  avg_eig_vector = as.matrix(avg_eig_vector)
  
  # get the current eigen vector corresponding to u(t)
  curr_eig_vector = principal_eig_vectors[,i]
  curr_eig_vector = as.matrix(curr_eig_vector)
  # calculate the dot-product of the avg and current eigen vector procured above
  z[, i] = 1 - crossprod(avg_eig_vector, curr_eig_vector)
}

# since we did not consider the first W' elements, we need to remove those corresponding 
# values from the z list above
z = z[-seq(1:w_dash)]

# threshold is calculated as median + (3 * moving_average)
z_median = median(z)
z_moving_average = calculate_moving_average(z)
threshold_value = z_median + (3 * z_moving_average)

# plot the points and anomalie
x11()
plot(z, xlab = "Day", ylab = "Z")
abline(h = threshold_value)

#Get z score index which are more than threshold values
anomalies = data.frame()
for(i in 1:length(z)) {
   if(z[i] >= threshold_value) {
    # adding w_dash since we removed the first W' entries previously
    # subtracting 1 because the day index starts from 0
    day = i + w_dash - 1
    anomalies = rbind(anomalies, c(day, z[i]))
   }
}
anomalies = anomalies[order(anomalies[,2],decreasing = TRUE),]

anomalies_count = dim(anomalies)[1]

# output_data contains the string that is to be written to a file as output
output_data = NA
if(anomalies_count <= 10) {
  output_data = c(anomalies_count, anomalies[,1]);
}else if(anomalies_count < 100) {
  output_data = c(anomalies_count, anomalies[1:10,1]);
} else {
  count = floor(anomalies_count/10)
  output_data = c(anomalies_count, anomalies[1:count,1]);
}

filename = paste("output", feature_type, sep = "_")
full_file_path = paste(output_folder_path, filename, sep = "/")

# open a file connection
file_conn = file(full_file_path)

output_data = as.character(output_data)

# write the content to the file
writeLines(output_data, file_conn)

# close the file connection
close(file_conn)

# snipper to persist the plots when run thro command line
message("Press Return To Continue")
invisible(readLines("stdin", n=1))

