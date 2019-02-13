# recoverIDs.R

recoverIDs <- function(ensp, mart = myMart) {
  # Purpose:
  #     Try to recover IDs for ensp to sym mapping from biomart
  # Parameters:
  #     ensp: (character) a vector of ensemble peptide IDs
  #     mart: (Mart)      an ensemble mart object
  #     "HGNC" must exist in the global namespace
  # Value:
  #     result: a dataframe with columns "ensp" containing the ensemble
  #             peptide IDs of the input that could be mapped, and "sym",
  #             which contains the corresponding HGNC symbols, and rownames
  #             ensp.

  # Note: to figure out the correct filters and attributes to use in
  #       querying a biomart, first fetch the filters and attributes with
  #       code like:
  #         filt <- biomaRt::listFilters(myMart)
  #         attr <- biomaRt::listAttributes(myMart)
  #       ... and then query:
  #         attr$name[grep("RefSeq", attr$description, ignore.case = TRUE)]

  # Define which attributes we want to fetch from biomart, and which columns
  # those match to in "HGNC":
  myAtt <- data.frame(biomart = c("uniprotswissprot",
                                  "refseq_mrna",
                                  "ucsc"),
                      HGNC =    c("UniProtID",
                                  "RefSeqID",
                                  "UCSCID"),
                      stringsAsFactors = FALSE)

  # Send off biomart query
  bm <- biomaRt::getBM(filters    =   "ensembl_peptide_id",
                       attributes = c("ensembl_peptide_id",
                                      myAtt$biomart),
                       values     = ensp,
                       mart       = myMart)

  if (nrow(bm) > 0) {                   # at least one match was returned
    bm$sym <- rep(NA, nrow(bm))         # add a column to hold map results
    for (iCol in seq_len(ncol(bm))) {   # replace all "" with NA
      # Careful: combining logical vectors that can include NA is tricky.
      # Select elements that are neither already NA nor not-empty
      sel <- ( ! is.na(bm[ , iCol])) & (bm[ , iCol] == "")
      bm[sel, iCol] <- NA               # replace
    }
  }

  for (iAtt in seq_len(nrow(myAtt))) { # iterate over all requested attributes
    thisBmAtt <- myAtt$biomart[iAtt]
    thisHuAtt <- myAtt$HGNC[iAtt]
    if ( ! all(is.na(bm[ , thisBmAtt]))) {
      # some IDs were returned
      IDs <- bm[ , thisBmAtt]
      # get the symbol for a match, NA otherwise
      syms <- HGNC$sym[match(IDs, HGNC[ , thisHuAtt], incomparables = NA)]
      sel <- ( ! is.na(syms))
      bm$sym[sel] <- syms[sel] # Overwrite those that are not NA.
      # If there are multiple IDs returned for one row
      # in effect we return the last one that was
      # matched.
    }
  }
  # Post-process. Careful: ensemble_peptide_ids are not necessarily
  # unique in biomart output.
  bm <- bm[! is.na(bm$sym), c("ensembl_peptide_id", "sym")]
  bm <- bm[! duplicated(bm$ensembl_peptide_id), ]
  matchedIDs <- match(ensp, bm$ensembl_peptide_id)

  esMap <- data.frame(ensp = ensp,
                      sym = bm$sym[matchedIDs],
                      stringsAsFactors = FALSE)
  rownames(esMap) <- esMap$ensp

  # drop NAs
  esMap <- esMap[ ! is.na(esMap$sym), ]

  return(esMap)
}


# ====  TESTS  =================================================================
if (FALSE) {
  # Enter your function tests here...

}


# [END]
