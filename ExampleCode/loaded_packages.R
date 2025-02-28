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
write.csv(combined_table, file = "~/package_list.csv", row.names = FALSE)



options(
     usethis.description = list(
          "Authors@R" = utils::person(
               "Meelad", "Amouzgar",
               email = "amouzgar@stanford.edu",
          ),

          Title = 'Hafez',
          Description = 'Landmark trajectory inference and time-series analysis.',
          Language =  "es"     )
)
usethis::use_description(fields = list(), check_name = TRUE, roxygen = TRUE)
usethis::use_gpl3_license()

