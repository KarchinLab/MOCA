

hierarchical_cluster_distance_matrix <- function(Matrix,
                                    ClusterMethod="average",
                                    DistanceMethod="euclidean",
                                    Dimension="Rows"
                                    )
    {

        if (Dimension == "Both") {
            row.order <- hclust(dist(Matrix, method=DistanceMethod), method=ClusterMethod)$order
            col.order <- hclust(dist(t(Matrix, method=DistanceMethod)), method=ClusterMethod)$order
            return(Matrix[row.order, col.order])
        }
        else if (Dimension == "Rows") {
            row.order <- hclust(dist(Matrix, method=DistanceMethod), method=ClusterMethod)$order
            return(Matrix[row.order,])
        }
        else if (Dimension == "Columns") {
            col.order <- hclust(dist(t(Matrix, method=DistanceMethod)), method=ClusterMethod)$order
            return(Matrix[, col.order])
        }
        else {
            stop("'Dimension' argument in 'hierarchical_cluster_distance_matrix'",
                 "must be set to either 'Columns', 'Rows', or 'Both'.")
        }
    }
