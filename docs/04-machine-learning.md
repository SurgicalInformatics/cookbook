# Machine learning

## Deep learning

### Pulling images from REDCap directly to argodeep

#### Original file names


```r
library(REDCapR)
uri = "https://redcap.cir.ed.ac.uk/api/"
token = "" # API token here
record_list = 1:318
field_list = c("photo", "photo_2", "photo_3", "photo_4")
event_list = c("wound_concerns_arm_2", "questionnaire_1_arm_2",
               "questionnaire_2_arm_2", "questionnaire_3_arm_2")
directory = "wound_raw" # destination directory must exist already


for(record in record_list){
  for(field in field_list){
    for(event in event_list){
      result = 
        tryCatch({      # suppress breaking error when no image in slot
          redcap_download_file_oneshot(
            record        = record,
            field         = field,
            redcap_uri    = uri,
            token         = token,
            event         = event,
            overwrite     = TRUE,
            directory     = directory
          )
        }, error=function(e){})
    }
  }
}
```

#### Named from REDCap record ID and event


```r
library(REDCapR)
uri = "https://redcap.cir.ed.ac.uk/api/"
token = "" # API token here
record_list = 1:318
field_list = c("photo", "photo_2", "photo_3", "photo_4")
event_list = c("wound_concerns_arm_2", "questionnaire_1_arm_2",
               "questionnaire_2_arm_2", "questionnaire_3_arm_2")
directory = "wound_named" # destination directory must exist already

for(record in record_list){
  for(field in field_list){
    for(event in event_list){
      file_name = paste0(record, "_", field, "_", event, ".jpg")
      result = 
        tryCatch({
          redcap_download_file_oneshot(
            record        = record,
            field         = field,
            redcap_uri    = uri,
            token         = token,
            event         = event,
            overwrite     = TRUE,
            directory     = directory,
            file_name     = file_name
          )
        }, error=function(e){})
    }
  }
}
```

