library(shiny)
library(survival)
library(DT)
library(rms)

### Read in pat file for options 
pat_template <- readRDS("data/fake_pat.rds")

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
          column(6, p("How would you like to intervene?")),
          column(6, selectInput("decision_intervention", NULL,
                                choices = c("Arbitrary intervention targeting all MFRs" = "adj_manual",
                                            "Arbitrary intervention targeting BMI" = "adj_bmi",
                                            "Arbitrary intervention targeting Non-HDL cholesterol" = "adj_nonhdl",
                                            "Arbitrary intervention targeting SBP" = "adj_sbp",
                                            "Smoking cessation" = "smoking_cessation",
                                            "Antihypertensive initiation" = "antihypertensives",
                                            "Statin initiation" = "statins")))
        ),
        uiOutput("conditional_input_decision_intervention")
      ),
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
        column(6, selectInput("ethnicity", NULL, choices = setNames(levels(pat_template$ethnicity), levels(pat_template$ethnicity))))
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
        column(6, p("H/O Rheumatoid arthritis:")),
        column(6, selectInput("ra", NULL, choices = c("Absent", "Present")))
      ),
      fluidRow(
        column(6, p("H/O systemic lupus erythematosus:")),
        column(6, selectInput("sle", NULL, choices = c("Absent", "Present")))
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
      DTOutput("riskTable"),
      textOutput("text_message")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  ### Conditional input
  output$conditional_input_decision_intervention <- renderUI({
    ### Arbitrary intervention targeting all MFRs
    if (input$decision_intervention == "adj_manual"){
      div(
        wellPanel(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("SBP at visit 0:")),
            column(6, numericInput("sbp", NULL, min = 30, max = 260, value = 140))
          ),
          fluidRow(
            column(6, p("Target SBP under intervention:")),
            column(6, numericInput("intervention_sbp", NULL, min = 30, max = 260, value = 120))
          )),
        wellPanel(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("BMI at visit 0:")),
            column(6, numericInput("bmi", NULL, min = 15, max = 50, value = 25))
          ),
          fluidRow(
            column(6, p("Target BMI under intervention:")),
            column(6, numericInput("intervention_bmi", NULL, min = 15, max = 50, value = 25))
          )),
        wellPanel(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("Non-HDL cholesterol at visit 0:")),
            column(6, numericInput("nonhdl", NULL, min = 0, max = 20, value = 3))
          ),
          fluidRow(
            column(6, p("Target Non-HDL cholesterol under intervention:")),
            column(6, numericInput("intervention_nonhdl", NULL,  min = 0, max = 20, value = 3))
          )),
        wellPanel(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("Smoking status at visit 0:")),
            column(6, selectInput("smoking", NULL, choices = c("Non-smoker", "Ex-smoker", "Current")))
          ),
          fluidRow(
            column(6, p("Target smoking status under intervention:")),
            column(6, selectInput("intervention_smoking", NULL, choices = c("Non-smoker", "Ex-smoker", "Current")))
          ))
      )
      ### Arbitrary intervention targeting BMI
    } else if (input$decision_intervention == "adj_bmi"){
      div(
        wellPanel(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("SBP at visit 0:")),
            column(6, numericInput("sbp", NULL, min = 30, max = 260, value = 140))
          )),
        wellPanel(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("BMI at visit 0:")),
            column(6, numericInput("bmi", NULL, min = 15, max = 50, value = 25))
          ),
          fluidRow(
            column(6, p("Target BMI under intervention:")),
            column(6, numericInput("intervention_bmi", NULL, min = 15, max = 50, value = 25))
          )),
        wellPanel(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("Non-HDL cholesterol at visit 0:")),
            column(6, numericInput("nonhdl", NULL, min = 0, max = 20, value = 3))
          )),
        wellPanel(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("Smoking status at visit 0:")),
            column(6, selectInput("smoking", NULL, choices = c("Non-smoker", "Ex-smoker", "Current")))
          ))
      )
      ### Arbitrary intervention targeting Non-HDL cholesterol
    } else if (input$decision_intervention == "adj_nonhdl"){
      div(
        wellPanel(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("SBP at visit 0:")),
            column(6, numericInput("sbp", NULL, min = 30, max = 260, value = 140))
          )),
        wellPanel(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("BMI at visit 0:")),
            column(6, numericInput("bmi", NULL, min = 15, max = 50, value = 25))
          )),
        wellPanel(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("Non-HDL cholesterol at visit 0:")),
            column(6, numericInput("nonhdl", NULL, min = 0, max = 20, value = 3))
          ),
          fluidRow(
            column(6, p("Target Non-HDL cholesterol under intervention:")),
            column(6, numericInput("intervention_nonhdl", NULL,  min = 0, max = 20, value = 3))
          )),
        wellPanel(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("Smoking status at visit 0:")),
            column(6, selectInput("smoking", NULL, choices = c("Non-smoker", "Ex-smoker", "Current")))
          ))
      )
      ### Arbitrary intervention targeting sbp
    } else if (input$decision_intervention == "adj_sbp"){
      div(
        wellPanel(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("SBP at visit 0:")),
            column(6, numericInput("sbp", NULL, min = 30, max = 260, value = 140))
          ),
          fluidRow(
            column(6, p("Target SBP under intervention:")),
            column(6, numericInput("intervention_sbp", NULL, min = 30, max = 260, value = 120))
          )),
        wellPanel(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("BMI at visit 0:")),
            column(6, numericInput("bmi", NULL, min = 15, max = 50, value = 25))
          )),
        wellPanel(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("Non-HDL cholesterol at visit 0:")),
            column(6, numericInput("nonhdl", NULL, min = 0, max = 20, value = 3))
          )),
        wellPanel(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("Smoking status at visit 0:")),
            column(6, selectInput("smoking", NULL, choices = c("Non-smoker", "Ex-smoker", "Current")))
          ))
      )
      ### Smoking cessation
    } else if (input$decision_intervention == "smoking_cessation"){
      wellPanel(
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("SBP at visit 0:")),
            column(6, numericInput("sbp", NULL, min = 30, max = 260, value = 140))
          )),
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("BMI at visit 0:")),
            column(6, numericInput("bmi", NULL, min = 15, max = 50, value = 25))
          )),
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("Non-HDL cholesterol at visit 0:")),
            column(6, numericInput("nonhdl", NULL, min = 0, max = 20, value = 3))
          )
        ),
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("Smoking status at visit 0:")),
            column(6, selectInput("smoking", NULL, choices = c("Current")))
          )
        )
      )
      ### Antihypertensives
    } else if (input$decision_intervention == "antihypertensives"){
      wellPanel(
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("SBP at visit 0:")),
            column(6, numericInput("sbp", NULL, min = 30, max = 260, value = 140))
          )),
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("BMI at visit 0:")),
            column(6, numericInput("bmi", NULL, min = 15, max = 50, value = 25))
          )),
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("Non-HDL cholesterol at visit 0:")),
            column(6, numericInput("nonhdl", NULL, min = 0, max = 20, value = 3))
          )
        ),
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("Smoking status at visit 0:")),
            column(6, selectInput("smoking", NULL, choices = c("Non-smoker", "Ex-smoker", "Current")))
          )
        )
      )
      ### Statins
    } else if (input$decision_intervention == "statins"){
      wellPanel(
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("SBP at visit 0:")),
            column(6, numericInput("sbp", NULL, min = 30, max = 260, value = 140))
          )),
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("BMI at visit 0:")),
            column(6, numericInput("bmi", NULL, min = 15, max = 50, value = 25))
          )),
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("Non-HDL cholesterol at visit 0:")),
            column(6, numericInput("nonhdl", NULL, min = 0, max = 20, value = 3))
          )
        ),
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("Smoking status at visit 0:")),
            column(6, selectInput("smoking", NULL, choices = c("Non-smoker", "Ex-smoker", "Current")))
          )
        )
      )
    }
  })
  
  ### Source functions
  source("data/functions_shiny.R")
  
  ### Read in pat file
  pat <- readRDS("data/fake_pat.rds")
  
  observeEvent(input$predictrisks, {
    
    if (input$gender == "Female"){
      fit_pui_list <- readRDS("data/fit_pui_female.rds")
      means_pui_list <- readRDS("data/means_pui_female.rds")
      bhaz_pui_list <- readRDS("data/bhaz_pui_female.rds")
      # fit_standard <- readRDS("data/fit_standard_female.rds")
      # bhaz_standard <- readRDS("data/bhaz_standard_female.rds")
    } else if (input$gender == "Male"){
      fit_pui_list <- readRDS("data/fit_pui_male.rds")
      means_pui_list <- readRDS("data/means_pui_male.rds")
      bhaz_pui_list <- readRDS("data/bhaz_pui_male.rds")
      # fit_standard <- readRDS("data/fit_standard_male.rds")
      # bhaz_standard <- readRDS("data/bhaz_standard_male.rds")
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
    pat$ra[1] <- input$ra
    pat$sle[1] <- input$sle
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
    
    ### Read in the direct effect hazard ratios
    direct_HR_smoking_initiation <- readRDS("data/offsets_direct_HR_smok_init.rds")
    direct_HR_smoking_cessation <- readRDS("data/offsets_direct_HR_smok_cess.rds")
    direct_HR_sbp <- readRDS("data/offsets_direct_HR_sbp.rds")
    direct_HR_bmi <- readRDS("data/offsets_direct_HR_bmi.rds")
    direct_HR_nonhdl <- readRDS("data/offsets_direct_HR_nonhdl.rds")
    
    ### Read in effect of mfrs on mfrs
    bmi_to_nonhdl <- readRDS("data/slope_bmi_to_nonhdl.rds")
    bmi_to_sbp <- readRDS("data/slope_bmi_to_sbp.rds")
    nonhdl_to_sbp <- readRDS("data/slope_nonhdl_to_sbp.rds")
    smok_cess_to_bmi <- readRDS("data/slope_smok_cess_to_bmi.rds")
    smok_cess_to_nonhdl <- readRDS("data/slope_smok_cess_to_nonhdl.rds")
    smok_cess_to_sbp <- readRDS("data/slope_smok_cess_to_sbp.rds")
    smok_init_to_bmi <- readRDS("data/slope_smok_init_to_bmi.rds")
    smok_init_to_nonhdl <- readRDS("data/slope_smok_init_to_nonhdl.rds")
    smok_init_to_sbp <- readRDS("data/slope_smok_init_to_sbp.rds")
    statins_to_bmi <- readRDS("data/slope_statins_to_bmi.rds")
    statins_to_nonhdl <- readRDS("data/slope_statins_to_nonhdl.rds")
    statins_to_sbp <- readRDS("data/slope_statins_to_sbp.rds")
    antihypertensives_to_sbp <- readRDS("data/slope_antihypertensives_to_sbp.rds")
    
    ###
    ### Estimate the relative risk change
    ###
    
    ### First get values at baseline
    baseline_smoking <- pat$smoking
    baseline_sbp <- pat$sbp
    baseline_bmi <- pat$bmi
    baseline_nonhdl <- pat$nonhdl
    
    ###  Arbitrary intervention targeting all risk factors
    if (input$decision_intervention == "adj_manual"){
      
      ### Assign values from intervention layer
      intervention_smoking <- input$intervention_smoking
      intervention_sbp <- input$intervention_sbp
      intervention_bmi <- input$intervention_bmi
      intervention_nonhdl <- input$intervention_nonhdl
      
      ### Apply floor values
      if (baseline_sbp < 120){
        baseline_sbp <- 120
      }
      if (intervention_sbp < 120){
        intervention_sbp <- 120
      }
      if (baseline_bmi < 25){
        baseline_bmi <- 25
      }
      if (intervention_bmi < 25){
        intervention_bmi <- 25
      }
      if (baseline_nonhdl < 2.6){
        baseline_nonhdl <- 2.6
      }
      if (intervention_nonhdl < 2.6){
        intervention_nonhdl <- 2.6
      }
      
      ### Get the difference
      change_sbp <- intervention_sbp - baseline_sbp
      change_bmi <- intervention_bmi - baseline_bmi
      change_nonhdl <- intervention_nonhdl - baseline_nonhdl
      
      ### Apply the relative risk change
      HR_change_sbp <- direct_HR_sbp^change_sbp
      HR_change_bmi <- direct_HR_bmi^change_bmi
      HR_change_nonhdl <- direct_HR_nonhdl^change_nonhdl
      
      ### Create output for text message
      output$text_message <- renderText({
        ""
      })
      
      ### Also need to do the same for smoking, which is different because it's non-continuous
      if (baseline_smoking == intervention_smoking){
        HR_change_smoking <- 1
      } else if (baseline_smoking == "Ex-smoker" & intervention_smoking == "Current"){
        HR_change_smoking <- 1/direct_HR_smoking_cessation
      } else if (baseline_smoking == "Current" & intervention_smoking == "Ex-smoker"){
        HR_change_smoking <- direct_HR_smoking_cessation
      } else if (baseline_smoking == "Non-smoker" & intervention_smoking == "Current"){
        HR_change_smoking <- direct_HR_smoking_initiation
      }  else if (baseline_smoking == "Non-smoker" & intervention_smoking == "Ex-smoker"){
        HR_change_smoking <- direct_HR_smoking_initiation*direct_HR_smoking_cessation
      } else if (baseline_smoking %in% c("Ex-smoker", "Current") & intervention_smoking == "Non-smoker"){
        HR_change_smoking <- NA
        ### Create output for text message
        output$text_message <- renderText({
          "You cannot move from a current or ex-smoker to non-smoker"
        })
      }
      # end undefined intervention block
    }
    
    ###  If intervention choice is to target BMI
    if (input$decision_intervention == "adj_bmi"){
      
      ### Assign values from intervention layer
      intervention_bmi <- input$intervention_bmi
      
      ### Calculate difference before applying floor values
      change_bmi_no_floor <- intervention_bmi - baseline_bmi
      
      ### Calculate expected knock-on changes before applying benefit floors
      change_nonhdl_no_floor <- bmi_to_nonhdl*change_bmi_no_floor
      change_sbp_no_floor <- bmi_to_sbp*change_bmi_no_floor
      
      ### Apply floor values
      if (baseline_bmi < 25){
        baseline_bmi <- 25
      }
      if (intervention_bmi < 25){
        intervention_bmi <- 25
      }
      
      ### Get the difference
      change_bmi <- intervention_bmi - baseline_bmi
      
      ### Calculate expected knock-on changes
      # Calculate expected change
      change_nonhdl <- bmi_to_nonhdl*change_bmi
      change_sbp <- bmi_to_sbp*change_bmi
      
      # Calculate
      intervention_nonhdl <- baseline_nonhdl + change_nonhdl
      intervention_sbp <- baseline_sbp + change_sbp
      
      ### Apply floor values
      if (baseline_sbp < 120){
        baseline_sbp <- 120
      }
      if (intervention_sbp < 120){
        intervention_sbp <- 120
      }
      if (baseline_nonhdl < 2.6){
        baseline_nonhdl <- 2.6
      }
      if (intervention_nonhdl < 2.6){
        intervention_nonhdl <- 2.6
      }
      
      ### Recalculate change after applying floors
      change_nonhdl <- intervention_nonhdl - baseline_nonhdl
      change_sbp <- intervention_sbp - baseline_sbp
      
      ### Apply the relative risk change
      HR_change_sbp <- direct_HR_sbp^change_sbp
      HR_change_bmi <- direct_HR_bmi^change_bmi
      HR_change_nonhdl <- direct_HR_nonhdl^change_nonhdl
      
      ### Set HR other MFRs, remain unchanged
      HR_change_smoking <- 1
      
      ### Create output for text message
      output$text_message <- renderText({
        paste0("A change in BMI of ", round(change_bmi_no_floor,2), " on average results in a change in Non-HDL cholesterol of ", round(change_nonhdl_no_floor, 2), 
               "mmol/L, and a change in SBP of ", round(change_sbp_no_floor,2), "mmHg. The hyopthetical risk is estimated based on this assumption. 
               These are the expected changes before applying benefit floors. 
               The hypothetical risk calculation applies benefit floors of BMI 25 kg/m², non-HDL cholesterol 2.6mmol/L and SBP 120 mmHg.")
      })
    }
    
    ###  If intervention choice is to target Non-HDL cholesterol
    if (input$decision_intervention == "adj_nonhdl"){
      
      ### Assign values from intervention layer
      intervention_nonhdl <- input$intervention_nonhdl
      
      ### Calculate difference before applying floor values
      change_nonhdl_no_floor <- intervention_nonhdl - baseline_nonhdl
      
      ### Calculate expected knock-on changes before applying benefit floors
      change_sbp_no_floor <- nonhdl_to_sbp*change_nonhdl_no_floor
      
      ### Apply floor values
      if (baseline_nonhdl < 2.6){
        baseline_nonhdl <- 2.6
      }
      if (intervention_nonhdl < 2.6){
        intervention_nonhdl <- 2.6
      }
      
      ### Get the difference
      change_nonhdl <- intervention_nonhdl - baseline_nonhdl
      
      ### Calculate expected knock-on changes
      # Calculate expected change
      change_sbp <- nonhdl_to_sbp*change_nonhdl
      
      # Calculate
      intervention_sbp <- baseline_sbp + change_sbp
      
      ### Apply floor values
      if (baseline_sbp < 120){
        baseline_sbp <- 120
      }
      if (intervention_sbp < 120){
        intervention_sbp <- 120
      }
      
      ### Recalculate change after applying floors
      change_sbp <- intervention_sbp - baseline_sbp
      
      ### Apply the relative risk change
      HR_change_sbp <- direct_HR_sbp^change_sbp
      HR_change_nonhdl <- direct_HR_nonhdl^change_nonhdl
      
      ### Set HR other MFRs, remain unchanged
      HR_change_smoking <- 1
      HR_change_bmi <- 1
      
      ### Create output for text message
      output$text_message <- renderText({
        paste0("A change in Non-HDL cholesterol of ", round(change_nonhdl_no_floor, 2), 
               " on average results in a change in SBP of ", round(change_sbp_no_floor, 2), 
               "mmHg. The hyopthetical risk is estimated based on this assumption.
               These are the expected changes before applying benefit floors. 
               The hypothetical risk calculation applies benefit floors of BMI 25 kg/m², non-HDL cholesterol 2.6mmol/L and SBP 120 mmHg.")
      })
    }
    
    ###  If intervention choice is to target SBP
    if (input$decision_intervention == "adj_sbp"){
      
      ### Assign values from intervention layer
      intervention_sbp <- input$intervention_sbp
      
      ### Apply floor values
      if (baseline_sbp < 120){
        baseline_sbp <- 120
      }
      if (intervention_sbp < 120){
        intervention_sbp <- 120
      }
      
      ### Get the difference
      change_sbp <- intervention_sbp - baseline_sbp
      
      ### Apply the relative risk change
      HR_change_sbp <- direct_HR_sbp^change_sbp
      
      ### Set HR other MFRs, remain unchanged
      HR_change_smoking <- 1
      HR_change_bmi <- 1
      HR_change_nonhdl <- 1
      
      ### Create output for text message
      output$text_message <- renderText({
        paste0("")
      })
    }
    
    ### If intervention choice is smoking cessation
    if (input$decision_intervention == "smoking_cessation"){
      
      ### Assign values for intervention layer
      ### Set smoking status to ex-smoker
      intervention_smoking <- "Ex-smoker"
      ## Apply effects on sbp, bmi and nonhdl
      intervention_sbp <- baseline_sbp + smok_cess_to_sbp
      intervention_bmi <- baseline_bmi + smok_cess_to_bmi
      intervention_nonhdl <- baseline_nonhdl + smok_cess_to_nonhdl
      
      ### Apply floor values
      if (baseline_sbp < 120){
        baseline_sbp <- 120
      }
      if (intervention_sbp < 120){
        intervention_sbp <- 120
      }
      if (baseline_bmi < 25){
        baseline_bmi <- 25
      }
      if (intervention_bmi < 25){
        intervention_bmi <- 25
      }
      if (baseline_nonhdl < 2.6){
        baseline_nonhdl <- 2.6
      }
      if (intervention_nonhdl < 2.6){
        intervention_nonhdl <- 2.6
      }
      
      ### Get the difference
      change_sbp <- intervention_sbp - baseline_sbp
      change_bmi <- intervention_bmi - baseline_bmi
      change_nonhdl <- intervention_nonhdl - baseline_nonhdl
      
      ### Apply the relative risk change
      HR_change_sbp <- direct_HR_sbp^change_sbp
      HR_change_bmi <- direct_HR_bmi^change_bmi
      HR_change_nonhdl <- direct_HR_nonhdl^change_nonhdl
      
      ### Also need to do the same for smoking, which is different because it's non-continuous
      if (baseline_smoking != "Current"){
        HR_change_smoking <- NA
        ### Create output for text message
        output$text_message <- renderText({
          "You cannot recieve benefit from smoking cessation if you are not a current smoker"
        })
      } else {
        HR_change_smoking <- direct_HR_smoking_cessation
        ### Create output for text message
        output$text_message <- renderText({
          paste0("Smoking cessation on average will also result in a change in BMI of ", round(smok_cess_to_bmi,2), ",
          a change in non-HDL cholesterol of ", round(smok_cess_to_nonhdl,2), "mmol/L, and a change in SBP of ", round(smok_cess_to_sbp,2), "mmHg.
          The hyopthetical risk is estimated based on this assumption. 
          These are the expected changes before applying benefit floors. 
          The hypothetical risk calculation applies benefit floors of BMI 25 kg/m², non-HDL cholesterol 2.6mmol/L and SBP 120 mmHg.")
        })
      }
      # end smoking cessation block
    }
    
    ### If intervention choice is antihypertensives
    if (input$decision_intervention == "antihypertensives"){
      
      ### Assign values for intervention layer
      intervention_smoking <- baseline_smoking
      intervention_bmi <- baseline_bmi
      intervention_nonhdl <- baseline_nonhdl
      ## Apply expected effects on sbp
      intervention_sbp <- baseline_sbp + antihypertensives_to_sbp
      
      ### Apply floor values
      if (baseline_sbp < 120){
        baseline_sbp <- 120
      }
      if (intervention_sbp < 120){
        intervention_sbp <- 120
      }
      
      ### Get the difference
      change_sbp <- intervention_sbp - baseline_sbp
      
      ### Apply the relative risk change
      HR_change_sbp <- direct_HR_sbp^change_sbp
      
      ### Set HR other MFRs, remain unchanged
      HR_change_smoking <- 1
      HR_change_bmi <- 1
      HR_change_nonhdl <- 1
      
      ### Create output for text message
      output$text_message <- renderText({
        paste0("Taking antihypertensives will on average result in a change in systolic blood pressure
          of ", round(antihypertensives_to_sbp, 2), ". However, this varies considerably
          for different individuals, and based on the strength of the antihypertensive. To see the risk
          if able to reduce systolic blood pressure by a targeted amount, use the `Arbitrary intervention on SBP` option.
          These are the expected changes before applying benefit floors. 
          The hypothetical risk calculation applies benefit floors of BMI 25 kg/m², non-HDL cholesterol 2.6mmol/L and SBP 120 mmHg.")
      })
      
      # end antihypertensives block
    }
    
    ### If intervention choice is statins
    if (input$decision_intervention == "statins"){
      
      ### Assign values for intervention layer
      ## Apply effects on sbp, bmi and nonhdl
      intervention_sbp <- baseline_sbp + statins_to_sbp
      intervention_bmi <- baseline_bmi + statins_to_bmi
      intervention_nonhdl <- baseline_nonhdl + statins_to_nonhdl
      
      ### Apply floor values
      if (baseline_sbp < 120){
        baseline_sbp <- 120
      }
      if (intervention_sbp < 120){
        intervention_sbp <- 120
      }
      if (baseline_bmi < 25){
        baseline_bmi <- 25
      }
      if (intervention_bmi < 25){
        intervention_bmi <- 25
      }
      if (baseline_nonhdl < 2.6){
        baseline_nonhdl <- 2.6
      }
      if (intervention_nonhdl < 2.6){
        intervention_nonhdl <- 2.6
      }
      
      ### Get the difference
      change_sbp <- intervention_sbp - baseline_sbp
      change_bmi <- intervention_bmi - baseline_bmi
      change_nonhdl <- intervention_nonhdl - baseline_nonhdl
      
      ### Apply the relative risk change
      HR_change_sbp <- direct_HR_sbp^change_sbp
      HR_change_bmi <- direct_HR_bmi^change_bmi
      HR_change_nonhdl <- direct_HR_nonhdl^change_nonhdl
      
      ### Set HR for smoking status, remains unchanged
      HR_change_smoking <- 1
      
      ### Set HR other MFRs, remain unchanged
      output$text_message <- renderText({
        paste0("Statin initiation on average will also result in an change in BMI of ", round(statins_to_bmi,2), ",
          a change in non-HDL cholesterol of ", round(statins_to_nonhdl,2), "mmol/L, and a change in SBP of ", round(statins_to_sbp,2), "mmHg.
          The hyopthetical risk is estimated based on this assumption.
          These are the expected changes before applying benefit floors. 
          The hypothetical risk calculation applies benefit floors of BMI 25 kg/m², non-HDL cholesterol 2.6mmol/L and SBP 120 mmHg.")
      })
      
      # end statins block
    }
    
    ### Apply these HRs
    HR_change_total <-
      HR_change_smoking*
      HR_change_sbp*
      HR_change_bmi*
      HR_change_nonhdl
    
    ### Create output
    risk_table <- create_table(pat = pat,
                               fit_pui_list = fit_pui_list,
                               bhaz_pui_list = bhaz_pui_list,
                               means_pui_list = means_pui_list,
                               HR_change_total = HR_change_total)
    
    output$riskTable <- renderDT({risk_table})
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)
