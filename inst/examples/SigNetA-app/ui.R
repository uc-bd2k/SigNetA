
library(shiny)
library(shinyjs)
library(DT)
library(shinydashboard)

dashboardPage(title="SigNetA",

  skin = "yellow",

  dashboardHeader(title = img(src='finre.jpg')),
  
  
  dashboardSidebar(
   
    sidebarMenu(
     # menuItem("Dashboard", tabName = "dashboard"),
      #menuItem("Raw data", tabName = "rawdata"),
      #h4("File Upload"),
      tags$br(),
      tags$br(),
      tags$br(),
      tags$br(),
    
    box(title="Input Data",solidHeader = TRUE,width=13,status="danger",collapsible = TRUE,collapsed = TRUE,
      fileInput('file1', 'Choose file to upload',
                accept = c(
                  'text/csv',
                  'text/comma-separated-values',
                  'text/tab-separated-values',
                  'text/plain',
                  '.csv',
                  '.tsv',
                  '.xls'
                )),
       
     # h4("Load Example"),
      actionButton("load_example","Load Example")),
    selectInput("algorithm",
                label="select an algorithm",choices=list("RWR"=1,"ModBioNet"=4,"SortNet"=5),                  
                selected=1),
    selectInput("PPI",
                label="select a network",choices=list("STRING"=1,"Interactome"=2),                  
                selected = 1),
   # tags$br(),
     # h4("Download Network"),
   box(
     title="Network Customization",solidHeader = TRUE,width=13,status="danger",collapsible = TRUE,collapsed = TRUE,
     box(
       title = "Layouts", status = "warning",
       collapsible = TRUE,
       collapsed=TRUE,
       solidHeader = TRUE,
       width = 13,
       selectInput("layout",
                   label="select a layout",choices=list("Fruchterman Reingold"=1,"Grid"=3,"Kamada Kawai"=4,"Sphere"=5,"GraphOpt"=6),                  
                   selected = 1)
     ),
   ##  box(
     ##  title = "Genes", status = "warning",
      ## collapsible = TRUE,
      ## collapsed=TRUE,
       ##solidHeader = TRUE,
       ##width = 13,
       ##textInput("nodeSelect", "Highlight Gene",value=c(), width = NULL, placeholder = NULL)
       # sliderInput("nodeSize", "Node Size:", 25, 100, 25)
       #selectInput("Focus", "Focus on node :",choices=list("BRCA1"="BRCA1","ATP"="ATP"),selected = "BRCA1")
   # # ),
     box(
       title = "Physics", status="warning",
       collapsible = TRUE,
       collapsed=TRUE,
       width = 13,
       
       solidHeader = TRUE,
       selectInput("phyactive", "Activate Physics :", choices=list("Yes"=TRUE,"No"=FALSE),selected=FALSE),
       selectInput("solve", "Physics :", choices=list("barnesHut"="barnesHut","repulsion"="repulsion","forceAtlas2Based"="forceAtlas2Based"),selected="barnesHut"),
       sliderInput("grav", "Gravitational Constant:", -3000, 0, 50),
       sliderInput("centralGrav", "Central Gravity:", 0, 10, 0.05)
     )
#      box(
#        title = "Pathways", status="warning",
#        collapsible = TRUE,
#        collapsed=TRUE,
#        width = 13,
#        
#        solidHeader = TRUE,
#        actionButton("goAn","Connect Pathways to Genes")
#      )
     ),
   
  ## box(title="Selections",solidHeader = TRUE,width=13,status="danger",collapsible = TRUE,collapsed = TRUE,
      ## box(
       ##  title = "Algorithms", status = "warning",
       ##  collapsible = TRUE,
        ## collapsed=TRUE,
        ## solidHeader = TRUE,
        ## width = 13,
         ##selectInput("algorithm",
                  ##   label="select an algorithm",choices=list("RWR"=1,"Bionet"=2,"dmGWAS"=3,"ModBioNet"=4,"SortNet"=5),                  
                   ##  selected=1),
      ## ),
      
      ## box(
       ##  title = "Protein Protein Interaction", status = "warning",
        ## collapsible = TRUE,
        ## collapsed=TRUE,
         ##solidHeader = TRUE,
        ## width = 13,
       ##  selectInput("PPI",
                  ##   label="select a network",choices=list("Interactome"=1,"STRING"=2),                  
                  ##   selected = 1),
      ## )
  ## ),
    box(
      title="Downloads",solidHeader = TRUE,width=13,status="danger",collapsible = TRUE,collapsed = TRUE,
      tags$br(),
      downloadButton('downloadNetworkData', 'Download Network Data'),
     tags$br(),
   
     # actionLink("saveImage", "Download as PNG"),
    tags$br(),
      #h4("Download Enrichment Results"),
      downloadButton('downloadEnrich', 'Download Enrichment Data')),
     
 
 
   # actionButton("goAn","Draw Pathways"),
    
      
      #tags$hr(),
      width = 3
    )
      
  ),
  
  
  
  dashboardBody(
    useShinyjs(),
    tags$head(tags$style(HTML('

       img{
                            /*width: 330px;
                             margin: 0% -7%;*/
                              /*width: 115%;*/
                                width:150px;
                               margin-top: 0%;
                              /*margin-left: -7%;*/
                                margin-left:-5%
       }
#maindivplot{
height:73vh !important
}
#plot{
height:73vh !important
}
       .main-header .logo{
            height:5%;
            /*width:330px*/
              /*width:18%*/
             width:230px
       }
        .main-sidebar{
          /*top:60px*/
        }
.navbar.navbar-static-top{
/*display:none !important;*/
margin-left:330px
background-color:#ecf0f5 !important

}
.skin-yellow .main-header .navbar {
background-color:#ecf0f5 !important
}

.content .row{
margin-top: -4%;
    position: absolute;
z-index: 10000;
    margin-left: 40px;
  
}
.content .col-sm-8{
margin-left:20px
}
.main-header{
background-color:#ecf0f5
}
.main-header .sidebar-toggle{
background-color:#222d32 !important; /*#f8b818*/
/*padding: 30px 30px*/
/*#margin-left: 1%;
   #margin-top: 10px;*/
padding: 10px;
border-radius: 4px;
position: absolute;
    margin-top: 73px;
margin-left:-10px


}
.skin-yellow .main-header .navbar .sidebar-toggle:hover{
background-color:#222d32 !important;
}

.main-header .sidebar-toggle:after{
//content:"Options"

}
#loadmessage {
                 position: fixed;
                 top: 400px;
                 left: 400px;
                 width: 50%;
                 padding: 8px 0px 8px 0px;
                 text-align: center;
                 font-weight: bold;
                 font-size: 100%;
                 color: #FFFFFF;
                 background-color: #337ab7;
                 border-radius:25px;
                 
                 z-index: 105;
}

.shiny-output-error { visibility: hidden; }
.shiny-output-error:before { visibility: hidden; }

.skin-yellow .sidebar a {color:black}
.sidebar-menu{
text-align:center
}
#load_example{
margin: auto;
    margin-top: -5%;
}

/*sidebar*/
.sidebar-menu .box-body{
padding:0px !important;
}
.sidebar-menu .box{
border-radius:0px !important;
background:none !important;
}
#shiny-modal{
z-index:10000
}
/*visNetwork*/
div.vis-configuration.vis-config-item{}
div.vis-configuration.vis-config-item.vis-config-s2{}
div.vis-configuration.vis-config-item.vis-config-s3{}
.box.box-solid.box-danger{border:#000000}
.box.box-solid.box-danger>.box-header{background:#000000}


                              '))),
    tags$head(tags$script(src="script.js")),
    fluidRow(
      
      tabBox(
        #title = "First tabBox",
        # The id lets us use input$tabset1 on the server to find the current tab
        id = "tabs", width=11,title="Network Analysis Details",
        tabPanel( "About",
                  includeHTML("index.html")
                  
        ),
        tabPanel("Network",
                 value="cytonet",
                 
                 h4("Network"),
                 
                 conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                  tags$div("Loading...",id="loadmessage")
                 ),
                 
                 textOutput(outputId="input_error"),
                 uiOutput("general_ui")
                 # rcytoscapejsOutput("plot", height="800px",width="900px")
                 
                 
                 
                 
                 
        ),
        tabPanel("Network Data",
                 value="datTable",
                 
                 conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                  tags$div("Loading...",id="loadmessage")
                 ),
                 textOutput(outputId="input_error1"),
                 
                 DT::dataTableOutput('contents')
                 
        ),
        
        tabPanel( "Gene Ontology(GO) Enrichment Analysis",
                  value="irich",
                  
                  
                  conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                   tags$div("Loading...",id="loadmessage")
                  ),
                  textOutput(outputId="input_error2"),
                # sliderInput("enrichFDR","Set FDR for enrichment",min=0.05,max=1,value=0.1),
                  textInput("enrichFDR","Set FDR for enrichment",0.05),
                 actionButton("goAn","Connect Pathways to Genes"),
                  DT::dataTableOutput('enrichmentAnalysis')
                  
                  
        )
      )

    
    
  )
  )
 
  
)
