library(parallel)
library(future)

# monitoring checklist
# each job should report time to complete a fold
# - the number of workers available
# - the system details
# - the average gpu time

monitor.report_available_workers <- function() {
  workers <- availableWorkers()
  cat(sprintf("#workders/#availableCores/#totalCores: %d/%d/%d, workers:\n", length(workers), availableCores(), detectCores()))
  workers
}

monitor.system_details <- function() {
  system("lscpu | egrep 'CPU(s)|per core|per socket'", intern = TRUE)
}

monitor.get_average_cpu_usage <- function() {
  system("ps -C R -o %cpu,%mem,drs,pid,cmd", intern = TRUE)
}
