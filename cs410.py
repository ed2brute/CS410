###############################################################################
# Pathology report search script.
# 
# Extracts diagnosis, tissue committe, and SNOMED codes from scanned pathology 
# reports received via RightFax by converting the PDF image to text using OCR.
# 
# Next the script indexes the scraped text and provides search capabilities
# via a GUI. 
# 
#
# Author: Jared Coberly
# Version 1.0
###############################################################################

import time
from datetime import datetime
import string 

import pytesseract                      # OCR utility
from pdf2image import convert_from_path # converts PDF to image
from spellchecker import SpellChecker   # identifies gibberish in the OCR conversion
import os

from whoosh import index
from whoosh.fields import Schema, TEXT, ID
from whoosh.qparser import QueryParser
#from whoosh.analysis import StemmingAnalyzer

import tkinter
from tkinter import scrolledtext
from tkhtmlview import HTMLScrolledText
import html

# the pathology report
class PathReport:
    acc_area = ''   #accession area: CSP, CCY
    acc_year = -1   #accession year: 2024, etc.
    acc_number = -1 #accession number: 1234
    tc_code = -1    #tissue committe code: 1, 7
    diagnosis = ""
    file = ""       #PDF file name
    page_number = ""       #page number in the PDF file

    def AccStr(self):
        return str(self.acc_area) + ": " + str(self.acc_year)[2:] + " " + str(self.acc_number)
    #end AccStr
#end PathReport

def ReadGYNCytoReports(pages, reports, filepath):
    #reports = list()
    accessions = set()

    for pg_num, page in enumerate(pages, start=1):
        page_text = pytesseract.image_to_string(page)
        page_words = page_text.split()

        # crude hack: use 'Results:' and 'Molecular' to bound the diagnosis text
        molecular_index = page_words.index("Molecular") if "Molecular" in page_words else -1
        results_index = page_words.index("Results:") if "Results:" in page_words else -1

        if results_index > 0 and molecular_index > 0:
            path_report = PathReport()   
        
            # save the accession information
            csp_index = page_words.index("CCY") if "CCY" in page_words else -1
            path_report.acc_area = "CCY"
            path_report.acc_year = 2000 + int(page_words[csp_index+1])
            path_report.acc_number = int(page_words[csp_index+2])
            path_report.page_number = str(pg_num)
            path_report.file = filepath
        
            # the reference lab is bad about having duplicate reports in the fax
            # if we've seen this report already in this fax, skip it
            if path_report.acc_number not in accessions:                
                
                path_report.diagnosis += "PAP and/or Molecular Results: "
                for i in range(results_index+1, molecular_index):                    
                    if page_words[i] == "HPV.mRNA":
                        path_report.diagnosis += "HPV mRNA "
                    elif page_words[i] == "po":                        
                        path_report.diagnosis = path_report.diagnosis[:-1] + ". "
                    else:
                        path_report.diagnosis += page_words[i].translate(str.maketrans('','',":&-\'\"‘“")) + " "

                reports.append(path_report)
                accessions.add(path_report.acc_number)

                #debug output
                print("Found CCY " + str(path_report.acc_year)[2:] + " " + str(path_report.acc_number) + ", reading report...")
            else:
                print("Found CCY " + str(path_report.acc_year)[2:] + " " + str(path_report.acc_number) + ", skipping (duplicate report)")
            #end if not duplicate                
        #end accession
    #end for pages    
    #print(len(reports))
#end ReadGYNCytoReports

def ReadDermReports(pages, reports, filepath):
    accessions = set()

    #extract all text from the fax
    for pg_num, page in enumerate(pages, start=1):
        page_text = pytesseract.image_to_string(page)
        page_words = page_text.split()
        
        # crude hack: use 'DIAGNOSIS:' and 'Electronic' to bound the diagnosis text
        diagnosis_index = page_words.index("DIAGNOSIS:") if "DIAGNOSIS:" in page_words else -1
        electronic_index = page_words.index("Electronic") if "Electronic" in page_words else -1

        if diagnosis_index > 0 and electronic_index > 0:
            path_report = PathReport()           

            # save the accession information
            csp_index = page_words.index("CSP") if "CSP" in page_words else -1
            path_report.acc_area = "CSP"
            path_report.acc_year = 2000 + int(page_words[csp_index+1])
            path_report.acc_number = int((page_words[csp_index+2]).translate(str.maketrans('','','string.punctuation'))) #delete any punctuation
            path_report.page_number = str(pg_num)
            path_report.file = filepath

            # the ref lab is bad about having duplicate reports in the fax
            # if we've seen this report already in this fax, skip it
            if path_report.acc_number not in accessions:                
                
                for i in range(diagnosis_index+1, electronic_index-4):
                    path_report.diagnosis += page_words[i] + " "

                reports.append(path_report)
                accessions.add(path_report.acc_number)

                #debug output
                print("Found CSP " + str(path_report.acc_year)[2:] + " " + str(path_report.acc_number) + ", reading report...")
            else:
                print("Found CSP " + str(path_report.acc_year)[2:] + " " + str(path_report.acc_number) + ", skipping (duplicate report)")
            #end if not duplicate                
        #end accession
    #end for page

    #print(len(reports))
#end ReadDermReports

def ReadSurgReports(pages, reports, filepath):
    # reports can mix non-gyn cytology with surgical pathology
    # we need to keep two lists of reports     
    surg_reports = list()
    cyto_reports = list()
    supplimental_reports = list()

    accessions = set()
    sup_accessions = set()

    #extract all text from the fax
    for pgnum, page in enumerate(pages, start=1):
        page_text = pytesseract.image_to_string(page)
        page_words = page_text.split()

        #print(page_words)

        # non-gyn cytology
        if page_text.find("CCY") > 0: 
            path_report = PathReport()
            print("NonGYN")
        
        # mismatch repair protien IHC
        elif page_text.find("MISMATCH REPAIR PROTEIN IMMUNOHISTOCHEMISTRY") > 0:             
            
            # crude hack: use 'Interpretation' and 'Electronically' to bound the diagnosis text
            begin_index = page_words.index("Interpretation") if "Interpretation" in page_words else -1
            end_index = page_words.index("Electronically") if "Electronically" in page_words else -1

            if begin_index > 0 and end_index > 0:
                path_report = PathReport()

                # save the accession information
                csp_index = page_words.index("(CSP") if "(CSP" in page_words else -1
                path_report.acc_area = "CSP"
                path_report.page_number = str(pgnum)
                path_report.file = filepath
                
                #strip hyphens if their are any
                if page_words[csp_index+1].find('-') > 0:
                    path_report.acc_year = 2000 + int(page_words[csp_index+1][:2])
                    path_report.acc_number = int(page_words[csp_index+1][3:-1])
                else:
                    path_report.acc_year = 2000 + int(page_words[csp_index+1])
                    path_report.acc_number = (page_words[csp_index+2])[:-1] #strip the trailing space and parenthesis

                if path_report.acc_number not in accessions:
                    for i in range(begin_index, end_index-1):
                        path_report.diagnosis += page_words[i] + " "

                    supplimental_reports.append(path_report)
                    accessions.add(path_report.acc_number)

                    #debug output
                    print("Found supplimental report to: CSP " + str(path_report.acc_year)[2:] + " " + str(path_report.acc_number) + ", reading report...")
                else:
                    print("Found supplimental report to: CSP " + str(path_report.acc_year)[2:] + " " + str(path_report.acc_number) + ", skipping (duplicate report)")
                #endif
            #endif
        # surg path
        else:
        
            # crude hack: use 'ACCESSION' and 'Pathologist' to bound the diagnosis text
            accession_index = page_words.index("ACCESSION#") if "ACCESSION#" in page_words else -1
            primary_pathologist = page_words.index("Pathologist:") if "Pathologist:" in page_words else -1

            if accession_index > 0 and primary_pathologist > 0:
                path_report = PathReport()           

                # save the accession information
                csp_index = page_words.index("(CSP") if "(CSP" in page_words else -1
                path_report.acc_area = "CSP"
                path_report.page_number = str(pgnum)
                path_report.file = filepath
                
                #strip hyphens if their are any
                if page_words[csp_index+1].find('-') > 0:
                    path_report.acc_year = 2000 + int(page_words[csp_index+1][:2])
                    path_report.acc_number = int(page_words[csp_index+1][3:-2])
                else:
                    path_report.acc_year = 2000 + int(page_words[csp_index+1])
                    path_report.acc_number = (page_words[csp_index+2])[:-2] #strip the trailing space and parenthesis

                # the ref lab is bad about having duplicate reports in the fax
                # if we've seen this report already in this fax, skip it
                if path_report.acc_number not in accessions:
                    # capture the TC code
                    if page_text.find("TCO1:") > 0 or page_text.find("TC code 1") > 0 or page_text.find("TC:O1") > 0: 
                        path_report.tc_code = 1
                    elif page_text.find("TCO7:") > 0 or page_text.find("TCO07 :") > 0 or page_text.find("TC code 7") > 0 or page_text.find("TC:O7") > 0:
                        path_report.tc_code = 7

                    # sometimes the OCR adds extra words, but they're non-sensical
                    # use spellcheck to remove them
                    spell = SpellChecker()
                    misspelled = spell.unknown(page_words[accession_index+15:accession_index+16])
                    if len(misspelled) > 0:
                        accession_index += 1

                    # make easy formatting corrections where we can
                    for i in range(accession_index+15, primary_pathologist-1):
                        
                        #common OCR problem, missed the space
                        if page_words[i] == "Nota":
                            path_report.diagnosis += "Not a "    
                        
                        #push the comment down one line for style
                        elif page_words[i] == "COMMENT:" or page_words[i] == "TCO1:" or page_words[i] == "TCO7:" or page_words[i] == "TCO07 :":
                            path_report.diagnosis += "" + page_words[i] + " "
                        
                        else:
                            path_report.diagnosis += page_words[i] + " "
                    
                    surg_reports.append(path_report)
                    accessions.add(path_report.acc_number)

                    #debug output
                    print("Found CSP " + str(path_report.acc_year)[2:] + " " + str(path_report.acc_number) + ", reading report... TC: " + str(path_report.tc_code))
                else:
                    print("Found CSP " + str(path_report.acc_year)[2:] + " " + str(path_report.acc_number) + ", skipping (duplicate report)")
                #end if not duplicate                
            #end accession
        #end if cyto vers surg path
    #end for page

    for rep in surg_reports:
        #num_reports += 1
        reports.append(rep)

    for rep in cyto_reports:
        #num_reports += 1
        reports.append(rep)

    for rep in supplimental_reports:
        #num_reports += 1
        reports.append(rep)

    #print(len(reports))
#end ReadSurgReports()

def SendToFile(reports):
    print(r"Writing to logfile: C:\Users\VHACMOCoberj\Desktop\Reports\log_" + str(datetime.now().date()) + ".txt")
    
    f = open(r"C:\Users\VHACMOCoberj\Desktop\Reports\log_" + str(datetime.now().date()) + ".txt", "a")

    for report in reports:
        f.write( report.AccStr() + " " + report.diagnosis + "\n")        
    #end for

    f.close()        
#end SendToFile
            
def enter_function(value):
    searchButton.invoke()
    print("enter")
#end enter_function

def DoSearch():
    print("Do Search: " + queryField.get())
    text_widget.delete("1.0", tkinter.END) #clear result field in GUI
    
    q = qp.parse(u"" + queryField.get())    
    try:
        searcher = ix.searcher()
        results = searcher.search(q)        

        html_str = "<html><body>"
        for r in results:
            print(r["acc"] + " - " + r["diagnosis"])
            html_str += "<p><a href='"+ r["filepath"] + "#page=" + r["pagenum"] + "'>" + r["acc"] + "</a> (page " + r["pagenum"] + ") :: " + r["diagnosis"] + "</p>"
        
        html_str += "</body></html>"
        text_widget.set_html(html.unescape(html_str))  #bug in set_html, it escapes special characters, making it impossible to use the '#' sign in the url  

    finally:
        searcher.close()        
        
        

###############################################################################
# main 
###############################################################################

reports = list()
start_time = datetime.now()

# set the path to the image and text processing engines and the PDF files
poppler_engine_path = r"C:/Users/VHACMOCoberJ/AppData/Local/Poppler/poppler-24.07.0/Library/bin"
pytesseract.pytesseract.tesseract_cmd = (r"C:/Users/VHACMOCoberJ/AppData/Local/Programs/Tesseract-OCR/tesseract.exe")
reports_dir = r"C:/Users/VHACMOCoberJ/Desktop/Reports"

# read the reports directory
dir_list = os.listdir(reports_dir)

# add an 's' to the end of files if there is more than one...just makes the program look smarter
plural = ''
if len(dir_list) > 1 :
    plural = 's'
print("Begin extraction: " + str(len(dir_list)) + " file" + plural)

# open the files in the reports directory and see what type of reports we have
for report_file in dir_list:
    try:
        pdfpages = convert_from_path(reports_dir + "/" + report_file, poppler_path=poppler_engine_path)
    
        page_text = pytesseract.image_to_string(pdfpages[0])
        page_words = page_text.split()
        
        # crude hack: find a word only found in a particular type of report
        # use that the identify the reports in this file
        if page_text.find("CYTOLOGY") > 0:
            print("GYN Cytology - " + report_file)
            ReadGYNCytoReports(pdfpages, reports, reports_dir + "/" + report_file)        
        elif page_text.find("Indianapolis,") > 0:
            print("Dermatopathology - " + report_file)
            ReadDermReports(pdfpages, reports, reports_dir + "/" + report_file)        
        elif page_text.find("Denver,") > 0 and page_text.find("AmeriPath") > 0:
            print("Surgical pathology, non-GYN cytology, and supplimental reports - " + report_file)
            ReadSurgReports(pdfpages, reports, reports_dir + "/" + report_file)        
        else:
            print("Non Anatomic Pathology file encountered: " + report_file)
    except:
        print("Skipping  '" + report_file + "' -- not a PDF or other error")
#end for

#SendToFile(reports)

end_time = datetime.now()
delta_time = (end_time - start_time).total_seconds()

print("End extraction: " + str(len(dir_list)) + " file" + plural)
print("Total reports entered: " + str(len(reports)))
print("Total extraction time: " + str(round(delta_time,2)) + " seconds (" + str(round(delta_time/60, 2)) + " minutes)")
SendToFile(reports) #debug output to plain text file

# build document schema and index the extracted text
print("Indexing...")
schema = Schema(acc=ID(stored=True), diagnosis=TEXT(stored=True), filepath=TEXT(stored=True), pagenum=TEXT(stored=True))
ix = index.create_in(reports_dir, schema)
writer = ix.writer()

for r in reports:
    a = r.AccStr()
    d = r.diagnosis
    pg = r.page_number
    fp = r.file
    writer.add_document(acc=u"" + a, diagnosis=u"" + d, filepath=u"" + fp, pagenum=u"" + pg)
writer.commit()

qp = QueryParser("diagnosis", schema=ix.schema)

print("Ready!")

# Simple python GUI
window = tkinter.Tk()
window.title("Pathology Report Search")
queryLabel = tkinter.Label(window, text="Query")
queryLabel.grid(row=0, column=0)

queryField = tkinter.Entry(window, width=60)
queryField.grid(row=0, column=1)

searchButton = tkinter.Button(window, text="Search", command=DoSearch)
searchButton.grid(row=0, column=2)
window.bind('<Return>', enter_function)

text_widget = HTMLScrolledText(window)
text_widget.grid(row=2, columnspan=3)

window.mainloop()

###############################################################################
