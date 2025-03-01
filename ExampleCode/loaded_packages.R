library(hafez)# Get session information
session_info <- sessionInfo()

# Extract attached packages
attached_packages <- session_info$otherPkgs

# Extract packages loaded via namespace (not attached)
namespace_packages <- loadedNamespaces()

# Create a data frame for attached packages
attached_table <- data.frame(
     Package = names(attached_packages),
     Version = sapply(attached_packages, function(pkg) pkg$Version),
     # Source = "Attached",
     stringsAsFactors = FALSE
)

package_versions =c()
for(i in namespace_packages){
     package_versions=c(package_versions, as.character(packageVersion(i)))
}
# Create a data frame for namespace packages
namespace_table <- data.frame(
     Package = namespace_packages,
     Version = package_versions,
     stringsAsFactors = FALSE
)

# Combine both tables
combined_table <- rbind(attached_table, namespace_table)

# Print the combined table
print(combined_table)

# Write the combined data frame to a CSV file
write.csv(combined_table, file = "package_list.csv", row.names = FALSE)



# options(
     # usethis.description = list(
     #      "Authors@R" = utils::person(
     #           "Meelad", "Amouzgar",
     #           email = "amouzgar@stanford.edu",
     #      ),
     #
     #      Title = 'Hafez',
     #      Description = 'Landmark trajectory inference and time-series analysis.',
     #      Language =  "es"     )
# )

# Use usethis::use_package() for each package
devtools::document()
usethis::use_description(fields =list(
     "Authors" = utils::person(
          "Meelad", "Amouzgar",
          email = "amouzgar@stanford.edu",
     ),

     Title = 'Hafez',
     Description = 'Landmark trajectory inference and time-series analysis.',
     Language =  "es"     ), check_name = TRUE, roxygen = TRUE)
usethis::use_gpl3_license()
usethis::use_pipe()
usethis::use_package("ElPiGraph.R")
usethis::use_package("dplyr")
usethis::use_package("dtwclust")
usethis::use_package("dynutils")
usethis::use_package("emmeans")
usethis::use_package("ggplot2")
usethis::use_package("igraph")
usethis::use_package("magrittr")
usethis::use_package("multcomp")
usethis::use_package("sqldf")
usethis::use_package("tidyr")
usethis::use_package("RANN")
usethis::use_package("mgcv")
usethis::use_package("scales")
usethis::use_package("tibble")
# load('data/cytof_example_data.rda')
usethis::use_data(example_cytof_data, overwrite = T)


devtools::check()


