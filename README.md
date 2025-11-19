PubR
Tools for Searching, Summarizing, and Visualizing PubMed Data in R

PubR is an educational R package designed to simplify the process of querying, extracting, processing, and visualizing bibliographic information from PubMed, based on a user-provided keyword.
The package interacts with the NCBI E-utilities API through the rentrez package and was developed primarily for learning and teaching purposes.

Features

PubR provides functions that automatically connect to PubMed using a search term supplied by the user and generate structured outputs and visualizations. Main features include:

• Word Cloud Generation

Creates a word cloud based on the most frequent words found in PubMed article titles retrieved using the user’s keyword.

• Yearly Publication Trend Plot

Plots the number of articles published per year for a given search keyword

• Most Frequent Authors

Identifies and lists the authors who appear most frequently in PubMed articles related to the user’s keyword.

• Most Common Journals

Provides a ranked list of journals that publish the highest number of articles on the searched topic.

Installation

Currently, the package is available on GitHub.

devtools::install_github("IcaroScalisse/PubR")

Dependencies

PubR uses the following packages:

rentrez

httr, jsonlite

tm, wordcloud, RColorBrewer

ggplot2

dplyr

All dependencies are installed automatically.

Main Functions
1. pubmed_wordcloud()

Creates a word cloud based on the most frequent words in PubMed article titles.

Example
pubmed_wordcloud("Artificial Intelligence", n_results = 400)

This function:

Retrieves PubMed IDs in safe batches

Fetches article titles

Cleans text (case, punctuation, numbers, stopwords)

Computes word frequencies

Generates a stylized word cloud

pubmed_year_publication()

Plots the yearly distribution of publications for a given keyword.

Example

pubmed_year_publication("artificial intelligence", year_range = c(2010, 2024))

pubmed_top_authors()

Retrieves PubMed summaries and identifies the most frequent authors associated with a keyword.

Example

pubmed_top_authors("CRISPR", top_n = 30)

top_journals()

Displays the journals that most frequently publish articles related to the specified keyword.

Example

top_journals("immunogenetics", top_n = 15)

get_top_words_by_year()

Retrieves the most common title words from articles published in a specific year.

Example

get_top_words_by_year(2020, top_n = 30)

License

This package is licensed under the MIT License, allowing free use, modification, and distribution.

Data & Stopwords

The package includes a helper data file:

data/stopwords.R — this file defines a character vector named stopwords (or similar) that contains words excluded from textual analyses (e.g., common words like “the”, “and”, or other domain-specific stopwords).

This vector is used by functions that process titles and generate word frequency tables and word clouds.

You can customize the list by editing data/stopwords.R before building/installing the package, or by supplying your own stopword vector in your R session (for development/testing with devtools::load_all() or via function arguments if you add that feature).
