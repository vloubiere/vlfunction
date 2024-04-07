# Define control group
#' Title
#'
#' @param group Data table containing the groups for which controls have to be selected
#' @param control Data table containing the potential controls
#' @param act.col.name Name of the column containing activity values (has to be the same in the two tables)
#'
#' @return Returns a subset of the control data.table corresponding with the act.col.name values compared to group
#' @export
#'
#' @examples
#' No exmaple yet
vl_select_act_matched <- function(group,
                                  control,
                                  act.col.name)
{
  ctl <- data.table::copy(control)
  sel <- rep(NA, nrow(group))
  while(anyNA(sel))
  {
    .c <- sapply(group[[act.col.name]][is.na(sel)], function(x) which.min(abs(x-ctl[[act.col.name]])))
    .c[duplicated(.c)] <- NA
    sel[is.na(sel)] <- .c
    ctl[(sel), (act.col.name):= NA]
    print(paste0(sum(!is.na(sel)), "/", nrow(group)))
  }
  return(control[(sel)])
}