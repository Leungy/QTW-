setwd("G:/My Drive/SMU/SMU/QTW/Week_3")

library(XML)
ubase = "http://www.cherryblossom.org/"
url = paste(ubase, "results/2012/2012cucb10m-f.htm", sep = "")
doc = htmlParse(url)

preNode = getNodeSet(doc, "//pre")

txt = xmlValue(preNode[[1]])

nchar(txt)

substr(txt, 1, 50)

substr(txt, nchar(txt) - 50, nchar(txt))

els = strsplit(txt, "\\r\\n")[[1]]

length(els)

els[1:3]

els[ length(els) ]

extractResTable =
  # Retrieve data from web site, find preformatted text,
  # return as a character vector.
  function(url)
  {
    doc = htmlParse(url)
    preNode = getNodeSet(doc, "//pre")
    txt = xmlValue(preNode[[1]])
    els = strsplit(txt, "\r\n")[[1]]
    return(els)
  }

ubase = "http://www.cherryblossom.org/"
urls = paste(ubase, "results/", 1999:2012, "/",
             1999:2012, "cucb10m-f.htm", sep = "")

womenTables = lapply(urls, extractResTable)

womenURLs =
  c("results/1999/cb99f.html", "results/2000/Cb003f.htm", "results/2001/oof_f.html",
    "results/2002/ooff.htm", "results/2003/CB03-F.HTM",
    "results/2004/women.htm", "results/2005/CB05-F.htm",
    "results/2006/women.htm", "results/2007/women.htm",
    "results/2008/women.htm", "results/2009/09cucb-F.htm",
    "results/2010/2010cucb10m-f.htm",
    "results/2011/2011cucb10m-f.htm",
    "results/2012/2012cucb10m-f.htm")

urls = paste(ubase, womenURLs, sep = "")
urls[1:3]

womenTables = lapply(urls, extractResTable)
names(womenTables) = 1999:2012

sapply(womenTables, length)

extractResTable =
  # Retrieve data from web site,
  # find the preformatted text,
  # and return as a character vector.
  function(url, year = 1999)
  {
    doc = htmlParse(url)
    if (year == 2000) {
      # Get preformatted text from 4th font element
      # The top file is ill formed so the <pre> search doesn't work.
      ff = getNodeSet(doc, "//font")
      txt = xmlValue(ff[[4]])
      els = strsplit(txt, "\r\n")[[1]]
    }
    #else if (year == 2009) {
      # Get preformatted text from <div class="Section1"> element
      # Each line of results is in a <pre> element
      #div1 = getNodeSet(doc, "//div[@class='Section1']")
      #pres = getNodeSet(div1[[1]], "//pre")
      #els = sapply(pres, xmlValue)
    #}
    else if (year == 1999) { # have to add this else if statement
      # Get preformatted text from <pre> elements
      pres = getNodeSet(doc, "//pre")
      txt = xmlValue(pres[[1]])
      els = strsplit(txt, "\n")[[1]]   
    } 
    else {
      # Get preformatted text from <pre> elements
      pres = getNodeSet(doc, "//pre")
      txt = xmlValue(pres[[1]])
      els = strsplit(txt, "\r\n")[[1]]   
    } 
  }

years = 1999:2012
womenTables = mapply(extractResTable, url = urls, year = years)

names(womenTables) = years
sapply(womenTables, length)

womenTables$'1999'[1:10]
womenTables[[2]][1:10] 

#### Save the outputs
save(womenTables, file = "CBwomenTextTables.rda")

###########################################################
###########################################################