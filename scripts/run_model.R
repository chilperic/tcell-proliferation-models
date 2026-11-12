# Run from repository root after installing the package.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript scripts/run_model.R config.json output-directory")
if (dir.exists(args[2L]) || file.exists(args[2L])) stop("Choose a new output directory")
config <- jsonlite::fromJSON(args[1L])
params <- unlist(config$parameters)
if (config$model == "cyton") {
  result <- tcellprolif::simulate_cyton(params, config$times, config$initial_cells, config$max_generation)
} else if (config$model == "smith_martin") {
  if (is.null(config$history) || is.null(config$legacy_survival)) stop("Specify history and legacy_survival explicitly")
  result <- tcellprolif::simulate_smith_martin(params, config$times, config$initial_cells, config$max_generation,
                                            history=config$history, legacy_survival=config$legacy_survival)
} else stop("Unknown model")
dir.create(args[2L], recursive=TRUE)
utils::write.csv(data.frame(time=attr(result,"times"),result,check.names=FALSE),file.path(args[2L],"trajectory.csv"),row.names=FALSE)
manifest <- list(schema_version="1.0", configuration=config, R=R.version.string,
                 platform=R.version$platform, packages=lapply(c("tcellprolif","deSolve","jsonlite","digest"),function(x) list(name=x,version=as.character(utils::packageVersion(x)))),
                 configuration_sha256=digest::digest(file=args[1L],algo="sha256"),
                 trajectory_sha256=digest::digest(file=file.path(args[2L],"trajectory.csv"),algo="sha256"),
                 method=attr(result,"method"), history=attr(result,"history"), legacy_survival=attr(result,"legacy_survival"))
jsonlite::write_json(manifest,file.path(args[2L],"run.json"),pretty=TRUE,auto_unbox=TRUE)
