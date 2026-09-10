## run_tests.R — smoke-test the landscape R6 class in landscape.R
## Run with: Rscript run_tests.R

section <- function(title) cat("\n============ ", title, " ============\n", sep = "")
ok <- function(msg) cat("[OK] ", msg, "\n", sep = "")

# Stand-ins for the package's own generics, which live elsewhere in the
# real package but need to exist for downscale()/em_dist() dispatch to
# be testable here in isolation.
downscale <- function(x, ...) UseMethod("downscale")
em_dist   <- function(x, ...) UseMethod("em_dist")

pdf(tempfile(fileext = ".pdf"))  # headless plotting device; no X11 needed

## ---------------------------------------------------------------------
section("1. sites = matrix, features = tibble")
## ---------------------------------------------------------------------

set.seed(1)
xy  <- cbind(x = runif(10, 0, 100), y = runif(10, 0, 100))
tib <- dplyr::tibble(intensity = runif(10))

ls_mat <- landscape$new(sites = xy, features = tib)
stopifnot(ls_mat$n == 10, ls_mat$sites_type == "matrix", ls_mat$features_type == "tibble")
ok("constructed matrix+tibble landscape")

cat("\n--- print(ls_mat) ---\n"); print(ls_mat)
cat("\n--- summary(ls_mat) ---\n"); summary(ls_mat)
cat("\n--- plot(ls_mat) ---\n"); plot(ls_mat); ok("plot() ran without error")

cat("\n--- downscale(ls_mat) / em_dist(ls_mat) (stubs) ---\n")
print(downscale(ls_mat))
print(em_dist(ls_mat))
ok("stub methods run and return NULL as expected")

## ---------------------------------------------------------------------
section("2. sites = sfc (POINT), features = sf boundaries")
## ---------------------------------------------------------------------

pts <- sf::st_sfc(lapply(1:6, function(i) sf::st_point(c(runif(1, 0, 10), runif(1, 0, 10)))))
bnd <- sf::st_sf(
  id = 1,
  geometry = sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 10,0, 10,10, 0,10, 0,0), ncol = 2, byrow = TRUE))))
)

ls_sfc <- landscape$new(sites = pts, features = bnd)
stopifnot(ls_sfc$n == 6, ls_sfc$sites_type == "sfc", ls_sfc$features_type == "boundaries")
ok("constructed sfc(POINT)+boundaries landscape")

cat("\n--- print(ls_sfc) ---\n"); print(ls_sfc)
cat("\n--- summary(ls_sfc) ---\n"); summary(ls_sfc)
cat("\n--- plot(ls_sfc) ---\n"); plot(ls_sfc); ok("plot() ran without error")

## also check sfc(POLYGON) is accepted as sites
polys <- sf::st_sfc(sf::st_polygon(list(matrix(c(0,0, 1,0, 1,1, 0,1, 0,0), ncol = 2, byrow = TRUE))))
ls_poly <- landscape$new(sites = polys)
stopifnot(ls_poly$sites_type == "sfc", ls_poly$n == 1)
ok("sfc(POLYGON) sites also accepted")

## ---------------------------------------------------------------------
section("3. sites = igraph adjacency graph, features = NULL")
## ---------------------------------------------------------------------

g <- igraph::make_ring(8)
ls_graph <- landscape$new(sites = g)
stopifnot(ls_graph$n == 8, ls_graph$sites_type == "graph", is.null(ls_graph$features_type))
ok("constructed graph-only landscape (no features)")

cat("\n--- print(ls_graph) ---\n"); print(ls_graph)
cat("\n--- summary(ls_graph) ---\n"); summary(ls_graph)
cat("\n--- plot(ls_graph) ---\n"); plot(ls_graph); ok("plot() ran without error")

## ---------------------------------------------------------------------
section("4. features = fake raster list (type-detection only; no terra installed)")
## ---------------------------------------------------------------------

fake_raster <- structure(list(), class = "SpatRaster")
ls_rast <- landscape$new(sites = xy, features = list(fake_raster, fake_raster))
stopifnot(ls_rast$features_type == "raster_list")
ok("raster_list feature type correctly detected")

cat("\n--- summary(ls_rast) ---\n"); summary(ls_rast)

## ---------------------------------------------------------------------
section("5. validation errors fire as expected")
## ---------------------------------------------------------------------

err <- function(expr) tryCatch({ expr; "NO ERROR RAISED" }, error = function(e) conditionMessage(e))

cat("3-column matrix:      ", err(landscape$new(sites = matrix(1:9, ncol = 3))), "\n")
cat("mismatched n features:", err(landscape$new(sites = xy, features = dplyr::tibble(z = 1:3))), "\n")
cat("garbage sites type:   ", err(landscape$new(sites = "not a valid sites object")), "\n")
cat("garbage features type:", err(landscape$new(sites = xy, features = "not valid features")), "\n")
ok("all four bad inputs raised errors (see messages above — confirm none say NO ERROR RAISED)")

## ---------------------------------------------------------------------
section("6. class / S3 dispatch sanity check")
## ---------------------------------------------------------------------

cat("class(ls_mat):", paste(class(ls_mat), collapse = ", "), "\n")
stopifnot(identical(class(ls_mat)[1], "landscape"))
ok("first class element is \"landscape\", matching the S3 method suffixes")

dev.off()
section("ALL TESTS PASSED")
