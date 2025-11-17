# Pacotes necessários
if (!requireNamespace("litsearchr", quietly = TRUE)) install.packages("litsearchr")
if (!requireNamespace("tm", quietly = TRUE)) install.packages("tm")
if (!requireNamespace("SnowballC", quietly = TRUE)) install.packages("SnowballC")
if (!requireNamespace("wordcloud", quietly = TRUE)) install.packages("wordcloud")
if (!requireNamespace("rentrez", quietly = TRUE)) install.packages("rentrez")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")

library(litsearchr)
library(tm)
library(SnowballC)
library(wordcloud)
library(rentrez)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(httr)
library(jsonlite)

#' Generate a PubMed Word Cloud
#'
#' Searches PubMed for titles containing the selected keyword, cleans the text,
#' processes the most frequent terms and generates a word cloud.
#'
#' @param keyword Character. Search term used in PubMed.
#' @param n_results Integer. Maximum number of results to retrieve.
#' @param pause Numeric. Delay (in seconds) between API requests to avoid rate limits.
#'
#' @return A dataframe of processed words and frequencies (invisible).
#' @export


pubmed_wordcloud <- function(keyword, n_results = 300, pause = 1) {
  message("🔍 Searching PubMed for keyword: ", keyword)

  # --- Batch search for IDs ---
  ids <- c()
  for (i in seq(0, n_results - 1, 300)) {
    batch_size <- min(300, n_results - i)
    search_res <- tryCatch({
      entrez_search(
        db = "pubmed",
        term = keyword,
        retmax = batch_size,
        retstart = i
      )
    }, error = function(e) {
      message("❌ Connection error during PubMed search.")
      return(NULL)
    })

    if (!is.null(search_res)) {
      ids <- c(ids, search_res$ids)
    }

    message("   ➜ IDs collected so far: ", length(ids))
    Sys.sleep(pause)
  }

  if (length(ids) == 0) {
    message("⚠️ No results found.")
    return(invisible(NULL))
  }

  # --- Batch fetch summaries ---
  message("📄 Fetching article titles...")
  titles <- c()
  for (j in seq(1, length(ids), 200)) {
    batch <- ids[j:min(j + 199, length(ids))]

    summaries <- tryCatch({
      entrez_summary(db = "pubmed", id = batch)
    }, error = function(e) {
      message("❌ Error retrieving summaries.")
      return(NULL)
    })

    batch_titles <- sapply(summaries, function(x) x$title)
    titles <- c(titles, batch_titles)

    Sys.sleep(pause)
  }

  titles <- titles[!is.na(titles)]

  if (length(titles) == 0) {
    message("⚠️ No titles found.")
    return(invisible(NULL))
  }

  # --- Text preprocessing ---
  message("⚙️ Cleaning and processing text...")
  corpus <- Corpus(VectorSource(titles))
  corpus <- tm_map(corpus, content_transformer(tolower))
  corpus <- tm_map(corpus, removePunctuation)
  corpus <- tm_map(corpus, removeNumbers)

  # 👉 Remove only your custom stopwords
  corpus <- tm_map(corpus, removeWords, stopwords)

  corpus <- tm_map(corpus, stripWhitespace)

  dtm <- TermDocumentMatrix(corpus)
  mtx <- as.matrix(dtm)
  freq <- sort(rowSums(mtx), decreasing = TRUE)
  df_words <- data.frame(word = names(freq), freq = freq)

  message("📊 Words processed: ", nrow(df_words))

  # --- Force keyword to appear emphasized ---
  keyword_lower <- tolower(keyword)
  if (!(keyword_lower %in% df_words$word)) {
    df_words <- rbind(df_words, data.frame(word = keyword_lower, freq = 0))
  }

  df_words$freq[df_words$word == keyword_lower] <- max(df_words$freq) * 2

  # --- Generate word cloud ---
  message("🖼️ Generating word cloud...")
  set.seed(122)
  par(mar = c(0, 0, 0, 0))
  suppressWarnings({
    wordcloud(
      words = df_words$word,
      freq = df_words$freq,
      min.freq = 2,
      max.words = 120,
      random.order = FALSE,
      rot.per = 0.1,
      scale = c(3.2, 0.7),
      colors = brewer.pal(8, "Dark2")
    )
  })

  message("✅ Word cloud successfully generated!")
  invisible(df_words)
}

#' Yearly Distribution of PubMed Publications
#'
#' Retrieves PubMed records for a specific keyword and plots the number of
#' publications per year, optionally filtered by a given year range.
#'
#' @param keyword Character. Search term.
#' @param n_results Integer. Maximum number of PubMed entries to fetch.
#' @param pause Numeric. Delay in seconds between batched requests.
#' @param year_range Numeric vector of length 2. Optional filter (start, end).
#'
#' @return A ggplot object showing publication counts per year.
#' @export

pubmed_year_publication <- function(
    keyword,
    n_results = 1000,
    pause = 0.5,
    year_range = NULL
) {
  message("🔍 Searching PubMed for keyword: ", keyword)

  # --- Batch search for IDs ---
  ids <- c()
  for (i in seq(0, n_results - 1, 300)) {
    batch <- tryCatch({
      entrez_search(
        db = "pubmed",
        term = keyword,
        retmax = min(300, n_results - i),
        retstart = i
      )
    }, error = function(e) {
      message("❌ Error connecting to PubMed.")
      return(NULL)
    })

    if (!is.null(batch)) {
      ids <- c(ids, batch$ids)
    }

    message("   ➜ IDs collected: ", length(ids))
    Sys.sleep(pause)
  }

  if (length(ids) == 0) {
    message("⚠️ No results found.")
    return(invisible(NULL))
  }

  # --- Batch fetch summaries ---
  summaries <- list()
  for (j in seq(1, length(ids), 200)) {
    batch_ids <- ids[j:min(j + 199, length(ids))]

    batch_summary <- tryCatch({
      entrez_summary(db = "pubmed", id = batch_ids)
    }, error = function(e) {
      message("❌ Error fetching summaries.")
      return(NULL)
    })

    summaries <- c(summaries, batch_summary)
    Sys.sleep(pause)
  }

  # --- Extract publication years ---
  pub_years <- sapply(summaries, function(x) x$pubdate)
  pub_years <- gsub("^([0-9]{4}).*", "\\1", pub_years)
  pub_years <- pub_years[grepl("^[0-9]{4}$", pub_years)]

  if (length(pub_years) == 0) {
    message("⚠️ No publication year information found.")
    return(invisible(NULL))
  }

  df_years <- as.data.frame(table(pub_years))
  names(df_years) <- c("Year", "Publications")
  df_years$Year <- as.numeric(as.character(df_years$Year))

  # --- Apply year range filter ---
  if (!is.null(year_range)) {

    if (length(year_range) != 2)
      stop("year_range must be a vector of length 2, e.g. c(2010, 2020).")

    if (year_range[1] > year_range[2])
      stop("year_range must be in the form c(start_year, end_year).")

    start_year <- year_range[1]
    end_year <- year_range[2]

    message(paste0("📅 Filtering years between ", start_year, " and ", end_year, "..."))

    df_years <- df_years[df_years$Year >= start_year & df_years$Year <= end_year, ]

    # If filtering removes all data
    if (nrow(df_years) == 0) {
      message("⚠️ No publications found inside the selected year range.")
      return(invisible(NULL))
    }

    # Fill missing years in the range
    full_range <- data.frame(Year = seq(start_year, end_year))
    df_years <- merge(full_range, df_years, by = "Year", all.x = TRUE)
    df_years$Publications[is.na(df_years$Publications)] <- 0
  }

  df_years <- df_years[order(df_years$Year), ]

  message("📊 Total publications considered: ", sum(df_years$Publications))

  # --- Plot ---
  ggplot(df_years, aes(x = Year, y = Publications)) +
    geom_col(fill = "#2E86AB", alpha = 0.8) +
    geom_line(color = "#1B4F72", linewidth = 0.8) +
    geom_point(color = "#1B4F72", size = 2) +
    geom_text(
      aes(label = ifelse(Publications > 0, Publications, "")),
      vjust = -0.4,
      size = 3
    ) +
    theme_minimal(base_size = 13) +
    labs(
      title = paste("Yearly distribution of publications on:", keyword),
      x = "Publication year",
      y = "Number of articles"
    ) +
    scale_x_continuous(
      breaks = df_years$Year  # ensures only the filtered years appear
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}


#' Identify the Most Frequent Authors for a PubMed Keyword
#'
#' Retrieves PubMed entries for a keyword, extracts authors from summaries,
#' counts occurrences, and plots the most frequent ones.
#'
#' @param keyword Character. Search term.
#' @param max_results Integer. Maximum number of results to retrieve.
#' @param batch_size Integer. Batch size for summary retrieval.
#' @param top_n Integer. Number of top authors to display.
#'
#' @return A dataframe containing authors and frequencies (invisible).
#' @export

pubmed_top_authors <- function(keyword, max_results = 1000, batch_size = 300, top_n = 30) {

  message("🔍 Searching PubMed for keyword: ", keyword)

  all_ids <- c()
  retstart <- 0
  web_hist <- NULL

  # ---- Batch search with retry + exponential backoff ----
  while (length(all_ids) < max_results) {

    fetch_limit <- min(batch_size, max_results - length(all_ids))
    retries <- 1
    max_retries <- 5

    repeat {
      search <- tryCatch({
        entrez_search(
          db = "pubmed",
          term = keyword,
          retmax = fetch_limit,
          retstart = retstart,
          use_history = TRUE
        )
      }, error = function(e) NULL)

      if (!is.null(search) && length(search$ids) > 0) break

      wait <- 2^(retries - 1)
      message("⏳ Request failed. Waiting ", wait, "s (attempt ", retries, ")")
      Sys.sleep(wait)

      retries <- retries + 1
      if (retries > max_retries) break
    }

    if (is.null(search) || length(search$ids) == 0) {
      message("⚠️ No more IDs returned.")
      break
    }

    all_ids <- c(all_ids, search$ids)
    retstart <- retstart + fetch_limit

    web_hist <- search$web_history   # <-- keep history from ANY successful batch

    if (length(all_ids) >= max_results) break
  }

  if (length(all_ids) == 0) {
    message("⚠️ No PubMed entries found.")
    return(invisible(NULL))
  }

  message("📄 Fetching summaries for ", length(all_ids), " articles via WebEnv...")

  # ---- Fetch summaries in batches (avoid 414 error) ----
  all_summaries <- list()
  fetched <- 0

  while (fetched < length(all_ids)) {
    fetch_limit <- min(300, length(all_ids) - fetched)

    summaries <- entrez_summary(
      db = "pubmed",
      web_history = web_hist,
      retstart = fetched,
      retmax = fetch_limit
    )

    all_summaries <- c(all_summaries, summaries)
    fetched <- fetched + fetch_limit
  }

  # ---- Extract authors ----
  extract_authors <- function(x) {
    aut <- x$authors
    if (is.null(aut)) return(character(0))

    if (is.list(aut)) {
      return(sapply(aut, function(a) {
        if (is.list(a) && "name" %in% names(a)) return(a$name)
        if (is.character(a)) return(a)
        return(NA)
      }))
    }

    if (is.character(aut)) return(aut)

    return(character(0))
  }

  authors <- unlist(lapply(all_summaries, extract_authors))

  # ---- Clean ----
  authors <- trimws(authors)
  authors <- authors[
    !is.na(authors) &
      authors != "" &
      !authors %in% c("Author", "AUTHORS", "Authors")
  ]

  if (length(authors) == 0) {
    message("⚠️ No valid authors found.")
    return(invisible(NULL))
  }

  # ---- Count ----
  df <- as.data.frame(table(authors))
  df <- df |>
    arrange(desc(Freq)) |>
    slice_head(n = top_n)

  # ---- Plot ----
  p <- ggplot(df, aes(x = reorder(authors, Freq), y = Freq)) +
    geom_col(fill = "#0077B6") +
    coord_flip() +
    theme_minimal(base_size = 13) +
    labs(
      title = paste("Top", top_n, "Most Frequent Authors for:", keyword),
      x = "Author",
      y = "Number of Publications"
    )

  print(p)
  invisible(df)
}

#' Most Common Journals for a PubMed Keyword
#'
#' Retrieves PubMed entries using WebEnv history (avoiding URL overflow)
#' and returns a bar plot of the journals that publish most on the topic.
#'
#' @param keyword Character. Search term.
#' @param max_results Integer. Maximum number of PubMed IDs to retrieve.
#' @param top_n Integer. Number of journals to display.
#' @param batch_size Integer. Number of summaries to fetch per batch.
#'
#' @return A ggplot bar chart showing the most frequent journals.
#' @export

top_journals <- function(keyword, max_results = 5000, top_n = 15, batch_size = 500) {
  message("🔍 Searching PubMed for keyword: ", keyword)

  # --- Search IDs + get WebEnv / QueryKey ---
  search_res <- entrez_search(
    db = "pubmed",
    term = keyword,
    retmax = max_results,
    usehistory = TRUE
  )

  if (length(search_res$ids) == 0) {
    stop("No results found for this keyword.")
  }

  message("📄 Fetching summaries for ", length(search_res$ids), " articles...")

  # --- Batch fetch summaries ---
  all_summaries <- list()
  num_batches <- ceiling(length(search_res$ids) / batch_size)

  for (i in seq_len(num_batches)) {
    start <- (i - 1) * batch_size + 1
    end <- min(i * batch_size, length(search_res$ids))

    message("  ⏳ Batch ", i, "/", num_batches, " (IDs ", start, " to ", end, ")")

    batch <- entrez_summary(
      db = "pubmed",
      web_history = search_res$web_history,
      retstart = start - 1,
      retmax = batch_size
    )

    all_summaries <- c(all_summaries, batch)
    Sys.sleep(0.34)  # Prevent rate-limit (3 requests/s)
  }

  # --- Extract journal names ---
  journals <- unlist(lapply(all_summaries, function(x) x$fulljournalname))

  df <- as.data.frame(table(journals))
  df <- df[order(df$Freq, decreasing = TRUE), ]
  df <- head(df, top_n)

  # --- Plot ---
  ggplot(df, aes(x = reorder(journals, Freq), y = Freq)) +
    geom_col() +
    coord_flip() +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Top", top_n, "Journals for:", keyword),
      x = "Journal",
      y = "Number of Publications"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

#' Retrieve the Most Frequent Words in PubMed Titles for a Given Year
#'
#' Searches PubMed articles from a specific year, extracts title words,
#' cleans them, removes stopwords and returns the most frequent ones.
#'
#' @param year Numeric. The publication year (YYYY).
#' @param n_results Integer. Maximum number of articles to retrieve.
#' @param top_n Integer. Number of most frequent words to display.
#' @param pause Numeric. Time (in seconds) between batch API requests.
#'
#' @return A ggplot bar chart, invisibly returns dataframe of words and freq.
#' @export

get_top_words_by_year <- function(year, n_results = 1000, top_n = 15, pause = 1) {

  message("🔍 Searching PubMed for articles from ", year, "...")

  # ===== 1. Retrieve IDs in batches =====
  ids <- c()
  for (i in seq(0, n_results - 1, 300)) {
    search_res <- entrez_search(
      db = "pubmed",
      term = paste0(year, "[DP]"),
      retmax = min(300, n_results - i),
      retstart = i
    )
    ids <- c(ids, search_res$ids)
    message("   ➜ IDs collected: ", length(ids))
    Sys.sleep(pause)
  }

  # ===== 2. Fetch summaries =====
  message("📄 Fetching article titles...")
  titles <- c()
  for (j in seq(1, length(ids), 200)) {
    batch_ids <- ids[j:min(j + 199, length(ids))]
    summaries <- entrez_summary(db = "pubmed", id = batch_ids)
    batch_titles <- sapply(summaries, function(x) x$title)
    titles <- c(titles, batch_titles)
    Sys.sleep(pause)
  }

  # ===== 3. Cleaning =====
  titles <- tolower(titles)
  titles <- gsub("[^a-zA-Z ]+", " ", titles)
  titles <- gsub("\\s+", " ", titles)

  words <- unlist(strsplit(titles, " "))
  words <- words[nchar(words) > 2]

  # ===== 4. Remove stopwords =====
  words <- words[!words %in% stopwords]

  # ===== 5. Count and select top words =====
  df <- as.data.frame(table(words)) %>%
    arrange(desc(Freq)) %>%
    head(top_n)

  # ===== 6. Plot =====
  ggplot(df, aes(x = reorder(words, Freq), y = Freq)) +
    geom_col(fill = "#1E88E5") +
    geom_text(aes(label = Freq), hjust = -0.2, size = 4) +
    coord_flip() +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Top", top_n, "Words in PubMed Titles -", year),
      x = "Word",
      y = "Number of Articles"
    ) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
    expand_limits(y = max(df$Freq) * 1.1)
}

