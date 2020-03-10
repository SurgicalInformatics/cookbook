
# Programming in rlang


## rlang

### What is rlang?

rlang is a low-level programming API for R which the tidyverse uses (meaning it speaks to R in as R like way as possible, rather than a 'high-level' - high level is more user orientated and interpretable). It  enables you to extend what the tidyverse can do and adapt it for your own uses. It's particularly good to use if you're doing lots of more 'programming' type R work, for example, building a package, making a complex shiny app or writing functions. It might also be handy if you're doing lots of big data manipulation and want to manipulate different datasets in the same way, for example.

In this chapter we'll discuss some uses of it.

### Dynamic calling of variables using dplyr

In this example, say we have a tibble of variables, but we want to apply dynamic changes to it (so we feed R a variable, that can change, either using another function like `purr::map` or in a ShinyApp). In this instance, specifying each variable and each different possible consequence using different logical functions would take forever and be very clunky. So we can use rlang to simply put a dynamic variable/object through the same function.

We'll use an example where we want to summarise by different outcomes in a dynamic way.


```r
library(tidyverse)
#Make data - here nonsense numbers on deaths per 100k from guns say

example_data = tibble(countries = c('UK', 'USA', 'Pakistan', 'Mexico', 'Ireland', 'Estonia'),
                      region = c('Europe', 'Americas', 'Asia', 'Americas', 'Europe', 'Europe'),
                      Death_from_guns = c(1, 200, 150, 450, 3, 3.5),
                      Death_from_smoking = c(100, 300, 140, 150, 120, 300))


#example function for summarising using dynamic variables and bang bangs
#note metric must be numeric
summarise_feature = function(df, col_var, ...){
  require(tidyverse)
  
  wurly_curly = function(.){#The wurly_curly function makes things nicer
  require(rlang)
  !!quo_name(enquo(.))}

  summary_nm_sum <- paste0("metric_", wurly_curly(col_var))#The new LHS variable must have a different name from that on the RHS

  df %>%
    group_by(...) %>% 
    summarise(!!summary_nm_sum := sum({{col_var}}))
}


#output new variable
example_data %>% 
  summarise_feature(Death_from_guns, region) 

#How to dynamically change names using mutate and rlang
#Bang bang variables into mutate/tidyverse(so they can be dynamically changed)

#Select metric
metric = 'metric_Death_from_guns' #example but can be any named column

metric_mutate = paste0('prefix_', metric) #some reason doesn't take these within the mutate function so has to be outside

example_data %>% 
  summarise_feature(Death_from_guns, region) %>% 
  mutate(!!metric_mutate := !!rlang::sym(metric) * 2)#apply arbitary multiplication by two
```
