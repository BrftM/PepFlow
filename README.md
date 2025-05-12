
# TCRFlow

TCRFlow is an R package designed to facilitate peptide-TCR screening workflows, leveraging the **pepitope** library by Michael Schubert. The package includes a Shiny-based graphical interface for data analysis and visualization.

---

## Start the Application

### **Setup and Requirements**

Before running the application for the first time, make sure to complete the setup as described in below.

### **Step 1: Download the Repository**
You can download the package from the GitHub repository:
[TCRFlow GitHub](https://github.com/BrftM/prototype_shiny_R/archive/refs/heads/main.zip)

---

## Requirements

### **Step 2: Install pepitope library**
Start with the installation of the **pepitope** library:
[pepitope Installation](https://mschubert.github.io/pepitope/index.html)

### **Step 3: Install Dependencies**
<details>
<summary><strong>Required packages: (Click to expand)</strong></summary>

```r
# Installation of devtools
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

### Choose the most suitable option for you
# Option 1: Installation from ZIP file
install.packages("C:/Users/username/Downloads/TCRFlow.zip", repos = NULL, type = "source")

# Option 2: Installation via RStudio GUI
## Tools -> Install Packages -> Install from: Package Archive File (.tar.gz, .zip) -> choose TCRFlow.zip -> Install

# Option 3: Installation from GitHub in the same way as pepitope
devtools::install_github("BrftM/TCRFlow", dependencies = TRUE)


# Load package
library(TCRFlow)

# Start app
runShinyApp()
```
</details>

<details>
<summary><strong>(Optional) Manual Development Setup: (Click to expand)</strong></summary>

#### **Working Directory (Development Only)**
You do **not** need to be in the package's working directory to start the app if the package is installed.  
The function `runShinyApp()` will locate the app regardless of your current working directory.  

However, if you are developing or testing the app directly from the source (without installation), ensure your working directory is set to the package root:
```r
setwd("path/to/TCRFlow")
```


```r
# Required packages for shiny (not automatically installed with pepitope)
install.packages(c("shiny", "shinyjs", "shinyFiles"))

# Install Bioconductor package BSgenome
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

# Run the Shiny app from the local directory (development only)
shiny::runApp("inst/shinyapp")
```
</details>

---

## Step 4: Running the Application

After completing the setup, you can start the Shiny app from **any directory** with:
```r
library(TCRFlow)
runShinyApp()
```

Enjoy exploring TCR-peptide interactions with TCRFlow!
