# Resolve explicit YAML configuration placeholders for R analysis scripts.
release_root <- function() normalizePath(Sys.getenv("RELEASE_ROOT", getwd()), mustWork = TRUE)

local_paths <- function(config_file = Sys.getenv("LOCAL_PATH_CONFIG", file.path(release_root(), "config", "local.yaml"))) {
    if (requireNamespace("yaml", quietly = TRUE)) {
        values <- yaml::read_yaml(config_file)
    } else {
        lines <- trimws(readLines(config_file, warn = FALSE))
        lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
        if (any(!grepl("^[A-Z_]+:[[:space:]]*.+$", lines))) stop("Use a flat YAML string mapping or install yaml")
        keys <- sub(":.*$", "", lines)
        if (anyDuplicated(keys)) stop("Duplicate local configuration key")
        scalars <- trimws(sub("^[^:]+:", "", lines))
        scalars <- vapply(scalars, function(x) {
            if ((startsWith(x, '"') && endsWith(x, '"')) || (startsWith(x, "'") && endsWith(x, "'"))) x <- substr(x, 2L, nchar(x) - 1L)
            x
        }, character(1), USE.NAMES = FALSE)
        values <- as.list(setNames(scalars, keys))
    }
    values
}

release_expand <- function(text) {
    matches <- regmatches(text, gregexpr("\\$\\{[A-Z_]+\\}", text))[[1]]
    for (placeholder in unique(matches)) {
        key <- substr(placeholder, 3L, nchar(placeholder) - 1L)
        value <- if (key == "REPO_ROOT") release_root() else Sys.getenv(key, "")
        if (!nzchar(value)) value <- local_paths()[[key]]
        if (is.null(value) || !is.character(value) || length(value) != 1L || !nzchar(value)) stop("Missing local configuration: ", key)
        text <- gsub(placeholder, value, text, fixed = TRUE)
    }
    text
}

local_path <- function(key, ...) file.path(path.expand(release_expand(paste0("${", key, "}"))), ...)
