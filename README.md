# CS410 Final Project

Hospitals often employ third-party pathology services to make medical diagnoses on tissue removed during surgery. These pathology laboratories use a diversity of software to generate reports and most frequently the hospitals receive these reports as bundled digital faxes of scanned paper reports (PDFs containing images of paper reports). A hospital may receive multiple faxes per day of completed pathology reports. Finding specific patient or diagnostic information in this collection is often a manual process of opening a PDF (representing a single fax and containing numerous reports) and browsing it for the desired information. By utilizing existing toolkits and pipelines, I intended to develop a system that will allow a user to search this collection. Major steps will include crawling the PDF documents, leveraging optical character recognition to convert the report images to text, generation of indexes, and creation of a search interface. The system should execute with a reasonable run-time and present relevant results. Effectiveness will be demonstrated using standard metrics.

## Files
* cs410.py :: the project python application
* cs410demo.mp4 :: a video showing the application in use (view raw file to launch video)
* cs410final.pdf :: documentation of the project, its implementation, usage, and performance
