library(shiny)
library(survival)
library(DT)
library(rms)

# Define UI for the application
ui <- fluidPage(
  
  titlePanel("CHARIOT Causal Risk Prediction Model"),
  
  sidebarLayout(
    sidebarPanel(
      width = 6,
      wellPanel(
        style = "border: 2px solid blue; padding: 15px;",
        actionButton("predictrisks", "Predict risk")),
      wellPanel(
        width = 6,
        style = "border: 2px solid purple; padding: 15px;",
        fluidRow(
          column(6, p("Would you like to suggest an intervention?")),
          column(6, selectInput("decision_intervention", NULL, 
                                choices = c("Adjust modifiable risk factors manually", 
                                            "Smoking cessation")))
        ),
        uiOutput("conditional_input_decision_intervention")
      ),
      wellPanel(
        style = "border: 2px solid green; padding: 15px;",
        fluidRow(
          column(6, p("SBP:")),
          column(6, numericInput("sbp", NULL, min = 30, max = 260, value = 140))
        ),fluidRow(
          column(6, p("Target SBP under intervention:")),
          column(6, numericInput("intervention.sbp", NULL, min = 30, max = 260, value = 120))
        )),
      wellPanel(
        style = "border: 2px solid green; padding: 15px;",
        fluidRow(
          column(6, p("Non-HDL cholesterol:")),
          column(6, numericInput("nonhdl", NULL, min = 0, max = 20, value = 3))
        ),
        fluidRow(
          column(6, p("Target Non-HDL cholesterol under intervention:")),
          column(6, numericInput("intervention.nonhdl", NULL,  min = 0, max = 20, value = 3))
        )),
      wellPanel(
        style = "border: 2px solid green; padding: 15px;",
        fluidRow(
          column(6, p("BMI:")),
          column(6, numericInput("bmi", NULL, min = 15, max = 50, value = 24))
        ),
        fluidRow(
          column(6, p("Target BMI under intervention:")),
          column(6, numericInput("intervention.bmi", NULL, min = 15, max = 50, value = 24))
        )),
      wellPanel(
        style = "border: 2px solid green; padding: 15px;",
        fluidRow(
          column(6, p("Smoking status:")),
          column(6, selectInput("smoking", NULL, choices = c("Non-smoker", "Ex-smoker", "Smoker")))
        ),
        fluidRow(
          column(6, p("Target smoking status under intervention:")),
          column(6, selectInput("intervention.smoking", NULL, choices = c("Non-smoker", "Ex-smoker", "Smoker")))
        )),
      fluidRow(
        column(6, p("Gender:")),
        column(6, selectInput("gender", NULL, choices = c("Female", "Male")))
      ),
      fluidRow(
        column(6, p("Age:")),
        column(6, numericInput("age", NULL, min = 18, max = 85, value = 50))
      ),
      fluidRow(
        column(6, p("Ethnicity:")),
        column(6, selectInput("ethnicity", NULL, choices = c("white","bangladeshi","black african","black caribbean",
                                                             "chinese","indian","other asian","other ethnic","pakistani")))
      ),      
      fluidRow(
        column(6, p("Diabetes status:")),
        column(6, selectInput("diabetes", NULL, choices = c("Absent", "Type1", "Type2")))
      ),
      fluidRow(
        column(6, p("H/O Atrial Fibrillation:")),
        column(6, selectInput("af", NULL, choices = c("Absent", "Present")))
      ),
      fluidRow(
        column(6, p("H/O Hypertension:")),
        column(6, selectInput("hypertension", NULL, choices = c("Absent", "Present")))
      ),
      fluidRow(
        column(6, p("H/O Chronic Kidney Disease:")),
        column(6, selectInput("ckd", NULL, choices = c("Absent", "Present")))
      ),
      fluidRow(
        column(6, p("H/O Chronic Obstructive Pulmonary Disease:")),
        column(6, selectInput("copd", NULL, choices = c("Absent", "Present")))
      ),
      fluidRow(
        column(6, p("H/O Family history CVD:")),
        column(6, selectInput("fhcvd", NULL, choices = c("Absent", "Present")))
      ),
      fluidRow(
        column(6, p("H/O Serious Mental Illness:")),
        column(6, selectInput("smi", NULL, choices = c("Absent", "Present")))
      ),
      fluidRow(
        column(6, p("H/O Migraine:")),
        column(6, selectInput("migraine", NULL, choices = c("Absent", "Present")))
      ),
      fluidRow(
        column(6, p("H/O Oral Cancer:")),
        column(6, selectInput("oral_cancer", NULL, choices = c("Absent", "Present")))
      ),
      fluidRow(
        column(6, p("H/O Brain Cancer:")),
        column(6, selectInput("brain_cancer", NULL, choices = c("Absent", "Present")))
      ),
      fluidRow(
        column(6, p("H/O Lung Cancer:")),
        column(6, selectInput("lung_cancer", NULL, choices = c("Absent", "Present")))
      ),
      fluidRow(
        column(6, p("H/O Blood Cancer:")),
        column(6, selectInput("blood_cancer", NULL, choices = c("Absent", "Present")))
      ),
      fluidRow(
        column(6, p("H/O Intellectual Disability:")),
        column(6, selectInput("int_dis", NULL, choices = c("Absent", "Present")))
      ),
      fluidRow(
        column(6, p("H/O Downs Syndrome:")),
        column(6, selectInput("downs", NULL, choices = c("Absent", "Present")))
      ),
      fluidRow(
        column(6, p("H/O Pre-Eclamspia:")),
        column(6, selectInput("pre_eclampsia", NULL, choices = c("Absent", "Present")))
      ),
      fluidRow(
        column(6, p("H/O Postnatal Depression:")),
        column(6, selectInput("postnatal_depression", NULL, choices = c("Absent", "Present")))
      ),
      fluidRow(
        column(6, p("H/O Impotence:")),
        column(6, selectInput("impotence", NULL, choices = c("Absent", "Present")))
      ),
      
      fluidRow(
        column(6, p("Current Corticosteroid Use:")),
        column(6, selectInput("cortico", NULL, choices = c("Absent", "Present")))
      ),
      fluidRow(
        column(6, p("Current Atypical Antipsychotic Use:")),
        column(6, selectInput("antipsy", NULL, choices = c("Absent", "Present")))
      ),
      fluidRow(
        column(6, p("IMD:")),
        column(6, numericInput("IMD", NULL, min = 1, max = 20, value = 10))
      )
    ),
    
    # Show the prediction results
    mainPanel(
      width = 6,
      h3("Prediction Result"),
      DTOutput("riskTable")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  ### Source functions
  source("data/functions_shiny.R")
  
  ### Read in pat file
  pat <- readRDS("data/fake.pat.rds")
  
  observeEvent(input$predictrisks, {
    
    if (input$gender == "Female"){
      fit.prototype3 <- readRDS("data/fit.prototype3.female.rds")
      means.prototype3 <- readRDS("data/means.prototype3.female.rds")
      fit.standard <- readRDS("data/fit.standard.female.rds")
      bhaz.prototype3 <- readRDS("data/bhaz.prototype3.female.rds")
      bhaz.standard <- readRDS("data/bhaz.standard.female.rds")
    } else if (input$gender == "Male"){
      fit.prototype3 <- readRDS("data/fit.prototype3.male.rds")
      means.prototype3 <- readRDS("data/means.prototype3.male.rds")
      fit.standard <- readRDS("data/fit.standard.male.rds")
      bhaz.prototype3 <- readRDS("data/bhaz.prototype3.male.rds")
      bhaz.standard <- readRDS("data/bhaz.standard.male.rds")
    }
    
    ### Create pat file
    pat$age <- input$age
    pat$ethnicity[1] <- input$ethnicity
    pat$bmi <- input$bmi
    pat$sbp <- input$sbp
    pat$nonhdl <- input$nonhdl
    pat$smoking[1] <- input$smoking
    pat$cortico[1] <- input$cortico
    pat$antipsy[1] <- input$antipsy
    pat$diabetes[1] <- input$diabetes
    pat$hypertension[1] <- input$hypertension
    pat$af[1] <- input$af
    pat$ckd[1] <- input$ckd
    pat$copd[1] <- input$copd
    pat$fhcvd[1] <- input$fhcvd
    pat$smi[1] <- input$smi
    pat$migraine[1] <- input$migraine
    pat$oral_cancer[1] <- input$oral_cancer
    pat$brain_cancer[1] <- input$brain_cancer
    pat$lung_cancer[1] <- input$lung_cancer
    pat$blood_cancer[1] <- input$blood_cancer
    pat$int_dis[1] <- input$int_dis
    pat$downs[1] <- input$downs
    pat$pre_eclampsia[1] <- input$pre_eclampsia
    pat$postnatal_depression[1] <- input$postnatal_depression
    pat$impotence[1] <- input$impotence
    pat$IMD[1] <- input$IMD
    
    ###
    ### Estimate the relative risk change
    ###
    
    ### First get values at baseline
    baseline.smoking <- pat$smoking
    baseline.sbp <- pat$sbp
    baseline.bmi <- pat$bmi
    baseline.nonhdl <- pat$nonhdl
    
    ### Assign values from intervention layer
    intervention.smoking <- input$intervention.smoking
    intervention.sbp <- input$intervention.sbp
    intervention.bmi <- input$intervention.bmi
    intervention.nonhdl <- input$intervention.nonhdl
    # intervention.smoking <- "Non-smoker"
    # intervention.sbp <- 130
    # intervention.bmi <- 25
    # intervention.nonhdl <- 4
    
    ### Apply limits
    if (baseline.sbp < 120){
      baseline.sbp <- 120
    }
    if (intervention.sbp < 120){
      intervention.sbp <- 120
    }
    if (baseline.bmi < 25){
      baseline.bmi <- 25
    }
    if (intervention.bmi < 25){
      intervention.bmi <- 25
    }
    if (baseline.nonhdl < 4){
      baseline.nonhdl <- 4
    }
    if (intervention.nonhdl < 4){
      intervention.nonhdl <- 4
    }
    
    ### Get the difference
    change.sbp <- intervention.sbp - baseline.sbp
    change.bmi <- intervention.bmi - baseline.bmi
    change.nonhdl <- intervention.nonhdl - baseline.nonhdl
    
    ### Read in the relative risks
    direct_RR_smoking_initiation <- readRDS("data/direct_RR_smoking_initiation")
    direct_RR_smoking_cessation <- readRDS("data/direct_RR_smoking_cessation")
    direct_RR_sbp <- readRDS("data/direct_RR_sbp")
    direct_RR_bmi <- readRDS("data/direct_RR_bmi")
    direct_RR_nonhdl <- readRDS("data/direct_RR_nonhdl")
    
    ### Apply the relative risk change
    rr.change.sbp <- direct_RR_sbp^change.sbp
    rr.change.bmi <- direct_RR_bmi^change.bmi
    rr.change.nonhdl <- direct_RR_nonhdl^change.nonhdl
    
    ### Also need to do the same for smoking, which is different because it's non-continuous
    if (baseline.smoking == intervention.smoking){
      rr.change.smoking <- 1
    } else if (baseline.smoking == "Ex-smoker" & intervention.smoking == "Smoker"){
      rr.change.smoking <- 1/direct_RR_smoking_cessation
    } else if (baseline.smoking == "Smoker" & intervention.smoking == "Ex-smoker"){
      rr.change.smoking <- direct_RR_smoking_cessation
    } else if (baseline.smoking == "Non-smoker" & intervention.smoking == "Smoker"){
      rr.change.smoking <- direct_RR_smoking_initiation
    }  else if (baseline.smoking == "Non-smoker" & intervention.smoking == "Ex-smoker"){
      rr.change.smoking <- direct_RR_smoking_initiation*direct_RR_smoking_cessation
    }
    
    ### Apply these RRs
    rr.change.total <-
      rr.change.smoking*
      rr.change.sbp*
      rr.change.bmi*
      rr.change.nonhdl
    
    ### Create output
    risk.table <- create_table(pat, 
                               fit.prototype3 = fit.prototype3, 
                               bhaz.prototype3 = bhaz.prototype3, 
                               means.statins = means.prototype3[["statins"]],
                               means.ah = means.prototype3[["ah"]],
                               means.smoking1 = means.prototype3[["smoking1"]],
                               means.smoking2 = means.prototype3[["smoking2"]],
                               fit.standard = fit.standard, 
                               bhaz.standard = bhaz.standard,
                               rr.change.total = rr.change.total,
                               rr.change.smoking,
                               rr.change.sbp,
                               rr.change.bmi,
                               rr.change.nonhdl)
    
    output$riskTable <- renderDT({risk.table})
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
