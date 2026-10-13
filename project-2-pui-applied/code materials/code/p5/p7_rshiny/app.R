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
                                choices = c("Undefined intervention" = "adj_manual",
                                            "Smoking cessation" = "smoking_cessation",
                                            "Antihypertensives" = "antihypertensives")))
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
    if (input$decision_intervention == "adj_manual"){
      wellPanel(
        div(style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
            fluidRow(
              column(6, p("SBP at visit 1:")),
              column(6, numericInput("sbp", NULL, min = 30, max = 260, value = 140))
            ),fluidRow(
              column(6, p("Target SBP under intervention:")),
              column(6, numericInput("intervention_sbp", NULL, min = 30, max = 260, value = 120))
            )),
        div(wellPanel(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("BMI at visit 1:")),
            column(6, numericInput("bmi", NULL, min = 15, max = 50, value = 24))
          ),
          fluidRow(
            column(6, p("Target BMI under intervention:")),
            column(6, numericInput("intervention_bmi", NULL, min = 15, max = 50, value = 24))
          ))),
        div(
          wellPanel(
            style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
            fluidRow(
              column(6, p("Non-HDL cholesterol at visit 1:")),
              column(6, numericInput("nonhdl", NULL, min = 0, max = 20, value = 3))
            ),
            fluidRow(
              column(6, p("Target Non-HDL cholesterol under intervention:")),
              column(6, numericInput("intervention_nonhdl", NULL,  min = 0, max = 20, value = 3))
            ))
        ),
        div(
          wellPanel(
            style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
            fluidRow(
              column(6, p("Smoking status at visit 1:")),
              column(6, selectInput("smoking", NULL, choices = c("Non-smoker", "Ex-smoker", "Current")))
            ),
            fluidRow(
              column(6, p("Target smoking status under intervention:")),
              column(6, selectInput("intervention_smoking", NULL, choices = c("Non-smoker", "Ex-smoker", "Current")))
            ))
        )
      )
    } else if (input$decision_intervention == "smoking_cessation"){
      wellPanel(
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("SBP at visit 1:")),
            column(6, numericInput("sbp", NULL, min = 30, max = 260, value = 140))
          )),
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("BMI at visit 1:")),
            column(6, numericInput("bmi", NULL, min = 15, max = 50, value = 24))
          )),
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("Non-HDL cholesterol at visit 1:")),
            column(6, numericInput("nonhdl", NULL, min = 0, max = 20, value = 3))
          )
        ),
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("Smoking status at visit 1:")),
            column(6, selectInput("smoking", NULL, choices = c("Current")))
          )
        )
      )
    } else if (input$decision_intervention == "antihypertensives"){
      wellPanel(
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("SBP at visit 1:")),
            column(6, numericInput("sbp", NULL, min = 30, max = 260, value = 140))
          )),
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("BMI at visit 1:")),
            column(6, numericInput("bmi", NULL, min = 15, max = 50, value = 24))
          )),
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("Non-HDL cholesterol at visit 1:")),
            column(6, numericInput("nonhdl", NULL, min = 0, max = 20, value = 3))
          )
        ),
        div(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("Smoking status at visit 1:")),
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

    ###
    ### Estimate the relative risk change
    ###

    ### First get values at baseline
    baseline_smoking <- pat$smoking
    baseline_sbp <- pat$sbp
    baseline_bmi <- pat$bmi
    baseline_nonhdl <- pat$nonhdl

    ###  If intervention choice is to change modifiable risk factors manually
    if (input$decision_intervention == "adj_manual"){

      ### Assign values from intervention layer
      intervention_smoking <- input$intervention_smoking
      intervention_sbp <- input$intervention_sbp
      intervention_bmi <- input$intervention_bmi
      intervention_nonhdl <- input$intervention_nonhdl

      ### Apply limits
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

    ### If intervention choice is smoking cessation
    if (input$decision_intervention == "smoking_cessation"){

      ### Assign values for intervention layer
      ### Set smoking status to ex-smoker
      intervention_smoking <- "Ex-smoker"
      ## Set sbp to SBP - 3.5
      intervention_sbp <- baseline_sbp - 3.5
      ## Set bmi to bmi + 0.66
      intervention_bmi <- baseline_bmi + 0.66
      intervention_nonhdl <- baseline_nonhdl + 0.066

      ### Apply limits
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
          "Smoking cessation on average will also result in an increase in BMI of 0.66,
          an increase in non-HDL cholesterol of 0.066  mmol/L, and a reduction in SBP of 3.5mmHg.
          We present the combined effect of all these changes."
        })
      }
      # end smoking cessation block
    }

    ### If intervention choice is smoking cessation
    if (input$decision_intervention == "antihypertensives"){

      ### Assign values for intervention layer
      ### Set smoking status to ex-smoker
      intervention_smoking <- baseline_smoking
      ## Set sbp to SBP - 5.1
      intervention_sbp <- baseline_sbp - 5.1
      ## Set bmi and non-hdl to have no change
      intervention_bmi <- baseline_bmi
      intervention_nonhdl <- baseline_nonhdl

      ### Apply limits
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
      HR_change_smoking <- 1

      ### Create output for text message
      output$text_message <- renderText({
        "Taking antihypertensives will on average result in a reduction of systolic blood pressure
          of 5.1, and no reduction in body mass index or Non-HDL cholesterol. To see the relative risk
          reduction if able to reduce systolic blood pressure by more than 5.1, use the `undefined intervention` option"
      })

      # end antihypertensives block
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
