extract_from_nc <- function(netcdf, coordinates, variable = NULL) {

    layer <- terra::rast(netcdf)

    results <- terra::extract(layer, coordinates, ID = F)

    results_final <- cbind(coordinates, results)

    return(as.data.frame(results_final))
}


get_nearby <- function(netcdf, variable, coordinates, mode = "queen", verbose = TRUE) {

    if (verbose) cat("Getting nearby valid cells for", nrow(coordinates), "records\n")

    if (mode == "queen") {
        fadj = 4
    } else {
        fadj = matrix(nrow = 5, ncol = 5)
        fadj[] <- c(rep(1, 12), 0, rep(1, 12))
    }

    layer <- terra::rast(netcdf)
    layer <- terra::subset(layer, variable)

    adj_cells <- terra::adjacent(layer, terra::cellFromXY(layer, coordinates),
                    directions = fadj)
    adj_cells <- as.data.frame(adj_cells)

    adj_cells$ID <- seq_len(nrow(adj_cells))
    adj_cells <- adj_cells |>
        tidyr::pivot_longer(1:(ncol(adj_cells)-1), names_to = "index", values_to = "cell") |>
        dplyr::select(-index)

    xy_adj <- terra::extract(layer[[1]], adj_cells$cell, xy = T)
    xy_adj <- cbind(xy_adj, adj_cells)

    adj_cells_result <- tapply(seq_len(nrow(xy_adj)), xy_adj$ID, function(idx){
        target <- xy_adj[idx,]
        values_c <- target[,3]
        id <- target$ID[1]
        if (all(is.na(values_c))) {
            to_return <- data.frame(new_lon = coordinates$decimalLongitude[id],
                                    new_lat = coordinates$decimalLatitude[id],
                                    value = NA)
        } else if (sum(!is.na(values_c)) == 1) {
            to_return <- target[!is.na(values_c),1:3]
            colnames(to_return) <- c("new_lon", "new_lat", "value")
        } else {
            focalcell <- coordinates[id,1:2]
            possible <- target[!is.na(values_c),1:3]
            true_x <- focalcell[1,1]
            true_y <- focalcell[1,2]
            diff <- abs(possible$y - true_y) + abs(possible$x - true_x)
            diff[diff >= 360] <- diff[diff >= 360] - 360
            to_return <- possible[order(diff)[1],]
            colnames(to_return) <- c("new_lon", "new_lat", "value")
        }
        return(to_return)
    })
    adj_cells_result <- do.call("rbind", adj_cells_result)

    coordinates_idx <- bind_cols(
        coordinates[,c("decimalLongitude", "decimalLatitude")],
        adj_cells_result[,c("value", "new_lon", "new_lat")]
    )

    return(coordinates_idx)
}
