
# PepFlow

PepFlow is an R package designed to facilitate peptide-TCR screening workflows, leveraging the **pepitope** library. Informations about the workflow behind the application is provided in the [pepitope library documentation](https://mschubert.github.io/pepitope/index.html). The package includes a Shiny-based graphical interface for data analysis and visualization.


## Start the Application

### **Step 1: Install PepFlow**

We require R>=4.5.0, as this includes important fixes on Bioconductor.
 
```r
# Installation of remotes
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Installation from GitHub in the same way as pepitope
remotes::install_github("BrftM/PepFlow", dependencies = TRUE)

```

### Step 2: Running PepFlow Application

After completing the setup, you can start the Shiny app from **any directory** with:
```r
PepFlow::runShinyApp()
```

Enjoy exploring TCR-peptide interactions with PepFlow!
