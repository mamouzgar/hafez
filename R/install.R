# R/install.R
install_dependencies <- function() {

     ## install gcc
     gcc_check <- system("gcc --version", intern = TRUE)

     # If gcc is not installed, install it using brew
     if (length(gcc_check) == 0 || grepl("command not found", gcc_check[1])) {
          system("brew install gcc")
     } else {
          cat("gcc is already installed:\n", gcc_check, "\n")
     }

     ## setup to install distutils
     dir_path <- "~/.R"
     file_path <- file.path(dir_path, "Makevars")
     # Create the directory if it does not exist
     if (!dir.exists(dir_path)) {
          dir.create(dir_path)
     }
     # Create the file if it does not exist
     if (!file.exists(file_path)) {
          file.create(file_path)
     }
     # Define the lines to be added
     lines_to_add <- c(
          "FC = /opt/homebrew/Cellar/gcc/11.3.0_2/bin/gfortran",
          "F77 = /opt/homebrew/Cellar/gcc/11.3.0_2/bin/gfortran",
          "FLIBS = -L/opt/homebrew/Cellar/gcc/11.3.0_2/lib/gcc/11"
     )

     # Append the lines to the Makevars file
     write(lines_to_add, file = file_path, append = TRUE)

     ## install distutils and elpigraph forked repos
     remotes::install_github("mamouzgar/distutils", build_vignettes = FALSE)
     remotes::install_github("mamouzgar/ElPiGraph.R", build_vignettes = FALSE)

     # # Install terminal packages (Linux example)
     # if (Sys.info()["sysname"] == "Linux") {
     #      system("sudo apt-get install -y package3")  # Replace with your terminal package
     # }
     #
     # # Add more OS checks if necessary (e.g., macOS, Windows)
     # if (Sys.info()["sysname"] == "macOS") {
     #      system("sudo apt-get install -y package3")  # Replace with your terminal package
     # }

}

# Call the function during package loading
.onLoad <- function(libname, pkgname) {
     suppressMessages(install_dependencies())
     packageStartupMessage("Welcome to MyPackage! Enjoy using it.")

}
