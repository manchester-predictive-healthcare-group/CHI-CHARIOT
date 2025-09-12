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
          column(6, p("How would you like to intervene?")),
          column(6, selectInput("decision_intervention", NULL,
                                choices = c("Undefined intervention" = "adj_manual",
                                            "Smoking cessation" = "smoking_cessation",
                                            "Antihypertensives" = "antihypertensives",
                                            "Statins" = "statins")))
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
              column(6, numericInput("intervention.sbp", NULL, min = 30, max = 260, value = 120))
            )),
        div(wellPanel(
          style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
          fluidRow(
            column(6, p("BMI at visit 1:")),
            column(6, numericInput("bmi", NULL, min = 15, max = 50, value = 24))
          ),
          fluidRow(
            column(6, p("Target BMI under intervention:")),
            column(6, numericInput("intervention.bmi", NULL, min = 15, max = 50, value = 24))
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
              column(6, numericInput("intervention.nonhdl", NULL,  min = 0, max = 20, value = 3))
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
              column(6, selectInput("intervention.smoking", NULL, choices = c("Non-smoker", "Ex-smoker", "Current")))
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
        div(style = "border: 2px solid green; padding: 15px; margin-bottom: 3px;",
            fluidRow(
              column(6, p("SBP at visit 1:")),
              column(6, numericInput("sbp", NULL, min = 30, max = 260, value = 140))
            ),fluidRow(
              column(6, p("Target SBP under intervention:")),
              column(6, numericInput("intervention.sbp", NULL, min = 30, max = 260, value = 120))
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
    } else if (input$decision_intervention == "statins"){
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

    ### Read in the direct effects
    direct_OR_smoking_initiation <- readRDS("data/direct_OR_smoking_initiation.rds")
    direct_OR_smoking_cessation <- readRDS("data/direct_OR_smoking_cessation.rds")
    direct_OR_sbp <- readRDS("data/direct_OR_sbp.rds")
    direct_OR_bmi <- readRDS("data/direct_OR_bmi.rds")
    direct_OR_nonhdl <- readRDS("data/direct_OR_nonhdl.rds")

    ###
    ### Estimate the relative risk change
    ###

    ### First get values at baseline
    baseline.smoking <- pat$smoking
    baseline.sbp <- pat$sbp
    baseline.bmi <- pat$bmi
    baseline.nonhdl <- pat$nonhdl

    ###  If intervention choice is to change modifiable risk factors manually
    if (input$decision_intervention == "adj_manual"){

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

      ### Apply the relative risk change
      or.change.sbp <- direct_OR_sbp^change.sbp
      or.change.bmi <- direct_OR_bmi^change.bmi
      or.change.nonhdl <- direct_OR_nonhdl^change.nonhdl

      ### Create output for text message
      output$text_message <- renderText({
        "When specifying an undefined intervention, we assume the effect of this intervention on each of the modifiable risk
        factors is known, and the 'knock-on' effects of changes to the modifiable risk factors are not propagated through the model.
        For example, if specifying a reduction in BMI under intervention, the resulting change in non-HDL cholesterol and
        SBP will not be automatically applied in the risk under intervention. This will be reviewed in future iterations."
      })

      ### Also need to do the same for smoking, which is different because it's non-continuous
      if (baseline.smoking == intervention.smoking){
        or.change.smoking <- 1
      } else if (baseline.smoking == "Ex-smoker" & intervention.smoking == "Current"){
        or.change.smoking <- 1/direct_OR_smoking_cessation
      } else if (baseline.smoking == "Current" & intervention.smoking == "Ex-smoker"){
        or.change.smoking <- direct_OR_smoking_cessation
      } else if (baseline.smoking == "Non-smoker" & intervention.smoking == "Current"){
        or.change.smoking <- direct_OR_smoking_initiation
      }  else if (baseline.smoking == "Non-smoker" & intervention.smoking == "Ex-smoker"){
        or.change.smoking <- direct_OR_smoking_initiation*direct_OR_smoking_cessation
      } else if (baseline.smoking %in% c("Ex-smoker", "Current") & intervention.smoking == "Non-smoker"){
        or.change.smoking <- NA
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
      intervention.smoking <- "Ex-smoker"
      ## Set sbp to SBP - 3.5
      intervention.sbp <- baseline.sbp - 3.5
      ## Set bmi to bmi + 0.66
      intervention.bmi <- baseline.bmi + 0.66
      ## Set non-HDL cholesterol to nonhdl + 0.068
      intervention.nonhdl <- baseline.nonhdl + 0.068

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

      ### Apply the relative risk change
      or.change.sbp <- direct_OR_sbp^change.sbp
      or.change.bmi <- direct_OR_bmi^change.bmi
      or.change.nonhdl <- direct_OR_nonhdl^change.nonhdl

      ### Also need to do the same for smoking, which is different because it's non-continuous
      if (baseline.smoking != "Current"){
        or.change.smoking <- NA
        ### Create output for text message
        output$text_message <- renderText({
          "You cannot recieve benefit from smoking cessation if you are not a current smoker"
        })
      } else {
        or.change.smoking <- direct_OR_smoking_cessation
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
    if (input$decision_intervention == "statins"){

      ### Assign values for intervention layer
      intervention.sbp <- baseline.sbp - 2.62
      intervention.bmi <- baseline.bmi + 0.33
      intervention.nonhdl <- baseline.nonhdl + -1.3268

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

      ### Apply the relative risk change
      or.change.sbp <- direct_OR_sbp^change.sbp
      or.change.bmi <- direct_OR_bmi^change.bmi
      or.change.nonhdl <- direct_OR_nonhdl^change.nonhdl

      ### Also need to do the same for smoking, which is different because it's non-continuous
      or.change.smoking <- 1

      ### Create output for text message
      output$text_message <- renderText({
        "Taking statins will on average result in a reduction of systolic blood pressure
          of 2.62, a small increase in BMI of 0.33, and a reduction in non-HDL cholesterol of 1.3368.
          To see the relative risk reduction if able to further reduce these modifiable risk factors, please
        use the `undefined intervention` option"
      })

      # end statin block
    }

    ### If intervention choice is antihypertensives
    if (input$decision_intervention == "antihypertensives"){

      ### Assign values from intervention layer
      intervention.sbp <- input$intervention.sbp

      ### Apply limits
      if (baseline.sbp < 120){
        baseline.sbp <- 120
      }
      if (intervention.sbp < 120){
        intervention.sbp <- 120
      }

      ### Get the difference
      change.sbp <- intervention.sbp - baseline.sbp

      ### Apply the relative risk change
      or.change.sbp <- direct_OR_sbp^change.sbp
      or.change.bmi <- 1
      or.change.nonhdl <- 1
      or.change.smoking <- 1

      ### Create output for text message
      output$text_message <- renderText({
        "The effect of antihypertensives is hard to predict, depending on whether
        you are initiating antihypertensives for the first time, or increasing antihypertensive intensity.
        We suggest estimating the potential risk under intervention depending on a targeted level of SBP
        under the antihypertensive therapy"
      })
      # end antihypertensives block
    }

    ### Apply these RRs
    or.change.total <-
      or.change.smoking*
      or.change.sbp*
      or.change.bmi*
      or.change.nonhdl

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
                               or.change.total = or.change.total,
                               or.change.smoking,
                               or.change.sbp,
                               or.change.bmi,
                               or.change.nonhdl)

    output$riskTable <- renderDT({risk.table})
  })

}

# Run the application
shinyApp(ui = ui, server = server)
