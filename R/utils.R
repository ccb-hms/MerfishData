.getResource <- function(recs, suffix)
{
    suffix <- paste0(suffix, "$")
    ind <- grep(suffix, recs$title)
    id <- recs$ah_id[ind]
    return(recs[[id]])
}
