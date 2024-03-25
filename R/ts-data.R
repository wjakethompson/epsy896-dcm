library(tidyverse)
library(fs)
library(here)

library(taylor)
library(spotifyr)
access_token <- get_spotify_access_token()

popularity <- tribble(
  ~album_name,                           ~album_uri,
  "Taylor Swift",                        "7mzrIsaAjnXihW3InKjlC3",
  "Fearless",                            "43OpbkiiIxJO8ktIB777Nn",
  "Fearless (Taylor's Version)",         "4hDok0OAJd57SGIT8xuWJH",
  "Speak Now",                           "5EpMjweRD573ASl7uNiHym",
  "Speak Now (Taylor's Version)",        "5AEDGbliTTfjOB8TSm1sxt",
  "Red",                                 "1KlU96Hw9nlvqpBPlSqcTV",
  "Red (Taylor's Version)",              "6kZ42qRrzov54LcAk4onW9",
  "1989",                                "34OkZVpuzBa9y40DCy0LPR",
  "1989 (Taylor's Version)",             "64LU4c1nfjz1t4VnGhagcg",
  "reputation",                          "6DEjYFkNZh67HP7R9PSZvv",
  "Lover",                               "1NAmidJlEaVgA3MpcPFYGq",
  "folklore",                            "1pzvBxYgT6OVwJLtHkrdQK",
  "evermore",                            "6AORtDjduMM3bupSWzbTSG",
  "Midnights",                           "1fnJ7k0bllNfL1kVdNVW1A"
) |> 
  mutate(popularity = map_int(album_uri,
                          function(.x) {
                            album <- get_album(.x)
                            album$popularity
                          }))


here("figure", "images", "taylor") |> 
  dir_ls() |> 
  as_tibble() |> 
  rename(img = value) |> 
  mutate(img = str_replace_all(img, "^.*(?=figure)", ""),
         era = path_file(img),
         era = path_ext_remove(era),
         era = str_replace(era, "[0-9]{2}-", ""),
         era = str_replace_all(era, "-", " "),
         era = case_when(era %in% c("reputation", "folklore", "evermore") ~
                           era,
                         str_detect(era, "tv") ~ str_replace(str_to_title(era),
                                                             "Tv", "(TV)"),
                         TRUE ~ str_to_title(era)),
         
         
         album = case_when(str_detect(era, "\\(TV\\)") ~
                             str_replace(era, "TV", "Taylor's Version"),
                           era == "Debut" ~ "Taylor Swift",
                           TRUE ~ era)) |> 
  relocate(img, .after = last_col()) |> 
  left_join(select(popularity, -album_uri),
            join_by(album == album_name)) |> 
  left_join(select(taylor_albums, -c(ep, album_release)),
            join_by(album == album_name)) |> 
  write_rds(here("data", "taylor.rds"))
