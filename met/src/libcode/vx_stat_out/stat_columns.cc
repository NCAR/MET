// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
// ** Copyright UCAR (c) 1992 - 2021
// ** University Corporation for Atmospheric Research (UCAR)
// ** National Center for Atmospheric Research (NCAR)
// ** Research Applications Lab (RAL)
// ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
// *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

////////////////////////////////////////////////////////////////////////

using namespace std;

#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <cstdio>
#include <cmath>

#include "stat_columns.h"

#include "vx_util.h"
#include "vx_log.h"

////////////////////////////////////////////////////////////////////////

void parse_row_col(const char *col_name, int &r, int &c) {
   const char *ptr = (const char *) 0;

   // Parse Fi_Oj strings
   r = atoi(++col_name);

   if((ptr = strrchr(col_name, '_')) != NULL) c = atoi(++ptr);
   else {
      mlog << Error << "\nparse_row_col() -> "
           << "unexpected column name specified: \""
           << col_name << "\"\n\n";
      throw(1);
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void open_txt_file(ofstream *&out, const char *file_name) {

   // Create and open the output file stream
   out = new ofstream;
   out->open(file_name);

   if(!(*out)) {
      mlog << Error << "\nopen_txt_file()-> "
           << "can't open the output file \"" << file_name
           << "\" for writing!\n\n";
      exit(1);
   }

   out->setf(ios::fixed);

   return;
}

////////////////////////////////////////////////////////////////////////

void close_txt_file(ofstream *&out, const char *file_name) {

   // List the file being closed
   mlog << Debug(1)
        << "Output file: " << file_name << "\n";

   // Close the output file
   out->close();
   delete out;
   out = (ofstream *) 0;

   return;
}

////////////////////////////////////////////////////////////////////////

ConcatString append_climo_bin(const ConcatString &mask_name,
                              int i_bin, int n_bin) {

   if(n_bin == 1) return(mask_name);

   // Append the climo CDF bin number.
   ConcatString cs;
   cs << mask_name << "_BIN_";
   if(i_bin == -1) cs << "MEAN";
   else            cs << i_bin+1;

   return(cs);
}

////////////////////////////////////////////////////////////////////////

void write_header_row(const char **cols, int n_cols, int hdr_flag,
                      AsciiTable &at, int r, int c) {
   int i;

   // Write the header column names if requested
   if(hdr_flag) {
      for(i=0; i<n_header_columns; i++)
         at.set_entry(r, i+c, (string)hdr_columns[i]);

      c += n_header_columns;
   }

   // Write the columns names specific to this line type
   for(i=0; i<n_cols; i++)
     at.set_entry(r, i+c, (string)cols[i]);

   return;
}

////////////////////////////////////////////////////////////////////////

void write_mctc_header_row(int hdr_flag, int n_cat, AsciiTable &at,
                           int r, int c) {
   int i, j, col;
   ConcatString cs;

   // Write the header column names if requested
   if(hdr_flag) {
      for(i=0; i<n_header_columns; i++)
         at.set_entry(r, i+c, (string)hdr_columns[i]);

      c += n_header_columns;
   }

   // Write the columns names specific to the MCTC line type
   at.set_entry(r, c+0, (string)mctc_columns[0]);
   at.set_entry(r, c+1, (string)mctc_columns[1]);

   // Write Fi_Oj for each cell of the NxN table
   for(i=0, col=c+2; i<n_cat; i++) {
      for(j=0; j<n_cat; j++) {
         cs.format("F%i_O%i", i+1, j+1);
         at.set_entry(r, col, cs); // Fi_Oj
         col++;
      }
   }

   at.set_entry(r, col, (string)mctc_columns[3]);
   col++;

   return;
}

////////////////////////////////////////////////////////////////////////

void write_pct_header_row(int hdr_flag, int n_thresh, AsciiTable &at,
                          int r, int c) {
   int i, col;
   ConcatString tmp_str;

   // Write the header column names if requested
   if(hdr_flag) {
      for(i=0; i<n_header_columns; i++)
         at.set_entry(r, i+c, (string)hdr_columns[i]);

      c += n_header_columns;
   }

   // Write the columns names specific to the PCT line type
   at.set_entry(r, c+0, (string)pct_columns[0]);
   at.set_entry(r, c+1, (string)pct_columns[1]);

   // Write THRESH_i, OY_i, ON_i for each row of the Nx2 table
   for(i=0, col=c+2; i<n_thresh-1; i++) {

      tmp_str.format("%s%i", pct_columns[2], i+1);
      at.set_entry(r, col, tmp_str); // Threshold
      col++;

      tmp_str.format("%s%i", pct_columns[3], i+1);
      at.set_entry(r, col, tmp_str); // Event Count (OY)
      col++;

      tmp_str.format("%s%i", pct_columns[4], i+1);
      at.set_entry(r, col, tmp_str); // Non-Event Count (ON)
      col++;
   }

   // Write out the last threshold
   if(n_thresh >= 1) {
      tmp_str.format("%s%i", pct_columns[2], n_thresh);
      at.set_entry(r, col, tmp_str);    // Threshold
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_pstd_header_row(int hdr_flag, int n_thresh, AsciiTable &at,
                           int r, int c) {
   int i, col;
   char tmp_str[max_str_len];

   // Write the header column names if requested
   if(hdr_flag) {
      for(i=0; i<n_header_columns; i++)
         at.set_entry(r, i+c, (string)hdr_columns[i]);

      c += n_header_columns;
   }

   // Write the columns names specific to the PSTD line type
   at.set_entry(r, c+0,  (string)pstd_columns[0]);
   at.set_entry(r, c+1,  (string)pstd_columns[1]);
   at.set_entry(r, c+2,  (string)pstd_columns[2]);
   at.set_entry(r, c+3,  (string)pstd_columns[3]);
   at.set_entry(r, c+4,  (string)pstd_columns[4]);
   at.set_entry(r, c+5,  (string)pstd_columns[5]);
   at.set_entry(r, c+6,  (string)pstd_columns[6]);
   at.set_entry(r, c+7,  (string)pstd_columns[7]);
   at.set_entry(r, c+8,  (string)pstd_columns[8]);
   at.set_entry(r, c+9,  (string)pstd_columns[9]);
   at.set_entry(r, c+10, (string)pstd_columns[10]);
   at.set_entry(r, c+11, (string)pstd_columns[11]);
   at.set_entry(r, c+12, (string)pstd_columns[12]);
   at.set_entry(r, c+13, (string)pstd_columns[13]);
   at.set_entry(r, c+14, (string)pstd_columns[14]);
   at.set_entry(r, c+15, (string)pstd_columns[15]);
   at.set_entry(r, c+16, (string)pstd_columns[16]);

   // Write THRESH_i for each threshold
   for(i=0, col=c+17; i<n_thresh; i++) {

      snprintf(tmp_str, sizeof(tmp_str), "%s%i", pstd_columns[17], i+1);
      at.set_entry(r, col, (string)tmp_str); // Threshold
      col++;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_pjc_header_row(int hdr_flag, int n_thresh, AsciiTable &at,
                          int r, int c) {
   int i, col;
   char tmp_str[max_str_len];

   // Write the header column names if requested
   if(hdr_flag) {
      for(i=0; i<n_header_columns; i++)
         at.set_entry(r, i+c, (string)hdr_columns[i]);

      c += n_header_columns;
   }

   // Write the columns names specific to the PJC line type
   at.set_entry(r, c+0, (string)pjc_columns[0]);
   at.set_entry(r, c+1, (string)pjc_columns[1]);

   // Write THRESH_i, OY_TP_i, ON_TP_i, CALIBRATION_i, REFINEMENT_i,
   // LIKELIHOOD_i, and BASER_i for each row of the Nx2 table
   for(i=0, col=c+2; i<n_thresh-1; i++) {

      snprintf(tmp_str, sizeof(tmp_str), "%s%i", pjc_columns[2], i+1);
      at.set_entry(r, col, (string)tmp_str); // Threshold
      col++;

      snprintf(tmp_str, sizeof(tmp_str), "%s%i", pjc_columns[3], i+1);
      at.set_entry(r, col, (string)tmp_str); // Event Count/N (OY_TP)
      col++;

      snprintf(tmp_str, sizeof(tmp_str), "%s%i", pjc_columns[4], i+1);
      at.set_entry(r, col, (string)tmp_str); // Non-Event Count/N (ON_TP)
      col++;

      snprintf(tmp_str, sizeof(tmp_str), "%s%i", pjc_columns[5], i+1);
      at.set_entry(r, col, (string)tmp_str); // Calibration
      col++;

      snprintf(tmp_str, sizeof(tmp_str), "%s%i", pjc_columns[6], i+1);
      at.set_entry(r, col, (string)tmp_str); // Refinement
      col++;

      snprintf(tmp_str, sizeof(tmp_str), "%s%i", pjc_columns[7], i+1);
      at.set_entry(r, col, (string)tmp_str); // Likelihood
      col++;

      snprintf(tmp_str, sizeof(tmp_str), "%s%i", pjc_columns[8], i+1);
      at.set_entry(r, col, (string)tmp_str); // Base Rate
      col++;
   }

   // Write out the last threshold
   if(n_thresh >= 1) {
      snprintf(tmp_str, sizeof(tmp_str), "%s%i", pct_columns[2], n_thresh);
      at.set_entry(r, col, (string)tmp_str);    // Threshold
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_prc_header_row(int hdr_flag, int n_thresh, AsciiTable &at,
                          int r, int c) {
   int i, col;
   char tmp_str[max_str_len];

   // Write the header column names if requested
   if(hdr_flag) {
      for(i=0; i<n_header_columns; i++)
         at.set_entry(r, i+c, (string)hdr_columns[i]);

      c += n_header_columns;
   }

   // Write the columns names specific to the PRC line type
   at.set_entry(r, c+0, (string)prc_columns[0]);
   at.set_entry(r, c+1, (string)prc_columns[1]);

   // Write THRESH_i, PODY_i, POFD_i for each row of the Nx2 table
   for(i=0, col=c+2; i<n_thresh-1; i++) {

      snprintf(tmp_str, sizeof(tmp_str), "%s%i", prc_columns[2], i+1);
      at.set_entry(r, col, (string)tmp_str); // Threshold
      col++;

      snprintf(tmp_str, sizeof(tmp_str), "%s%i", prc_columns[3], i+1);
      at.set_entry(r, col, (string)tmp_str); // PODY
      col++;

      snprintf(tmp_str, sizeof(tmp_str), "%s%i", prc_columns[4], i+1);
      at.set_entry(r, col, (string)tmp_str); // POFD
      col++;
   }

   // Write out the last threshold
   if(n_thresh >= 1) {
      snprintf(tmp_str, sizeof(tmp_str), "%s%i", prc_columns[2], n_thresh);
      at.set_entry(r, col, (string)tmp_str);    // Threshold
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_eclv_header_row(int hdr_flag, int n_pnt, AsciiTable &at,
                           int r, int c) {
   int i, col;
   ConcatString tmp_str;

   // Write the header column names if requested
   if(hdr_flag) {
      for(i=0; i<n_header_columns; i++)
         at.set_entry(r, i+c, (string)hdr_columns[i]);

      c += n_header_columns;
   }

   // Write the columns names specific to the ECLV line type
   at.set_entry(r, c+0, (string)eclv_columns[0]);
   at.set_entry(r, c+1, (string)eclv_columns[1]);
   at.set_entry(r, c+2, (string)eclv_columns[2]);
   at.set_entry(r, c+3, (string)eclv_columns[3]);

   // Write CL_i and VALUE_i for each bin
   for(i=0, col=c+4; i<n_pnt; i++) {

      tmp_str.format("%s%i", eclv_columns[4], i+1);
      at.set_entry(r, col, tmp_str);
      col++;

      tmp_str.format("%s%i", eclv_columns[5], i+1);
      at.set_entry(r, col, tmp_str);
      col++;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_rhist_header_row(int hdr_flag, int n_rank, AsciiTable &at,
                            int r, int c) {
   int i, col;
   char tmp_str[max_str_len];

   // Write the header column names if requested
   if(hdr_flag) {
      for(i=0; i<n_header_columns; i++)
         at.set_entry(r, i+c, (string)hdr_columns[i]);

      c += n_header_columns;
   }

   // Write the columns names specific to the RHIST line type
   at.set_entry(r, c+0, (string)rhist_columns[0]);
   at.set_entry(r, c+1, (string)rhist_columns[1]);

   // Write RANK_i for each rank
   for(i=0, col=c+2; i<n_rank; i++) {

      snprintf(tmp_str, sizeof(tmp_str), "%s%i", rhist_columns[2], i+1);
      at.set_entry(r, col, (string)tmp_str); // Counts for each rank
      col++;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_phist_header_row(int hdr_flag, int n_bin, AsciiTable &at,
                            int r, int c) {
   int i, col;
   char tmp_str[max_str_len];

   // Write the header column names if requested
   if(hdr_flag) {
      for(i=0; i<n_header_columns; i++)
         at.set_entry(r, i+c, (string)hdr_columns[i]);

      c += n_header_columns;
   }

   // Write the columns names specific to the PHIST line type
   at.set_entry(r, c+0, (string)phist_columns[0]);
   at.set_entry(r, c+1, (string)phist_columns[1]);
   at.set_entry(r, c+2, (string)phist_columns[2]);

   // Write BIN_i for each bin
   for(i=0, col=c+3; i<n_bin; i++) {

      snprintf(tmp_str, sizeof(tmp_str), "%s%i", phist_columns[3], i+1);
      at.set_entry(r, col, (string)tmp_str); // Counts for each bin
      col++;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_orank_header_row(int hdr_flag, int n_ens, AsciiTable &at,
                            int r, int c) {
   int i, col;
   char tmp_str[max_str_len];

   // Write the header column names if requested
   if(hdr_flag) {
      for(i=0; i<n_header_columns; i++)
         at.set_entry(r, i+c, (string)hdr_columns[i]);

      c += n_header_columns;
   }

   // Write the columns names specific to the orank line type
   at.set_entry(r, c+0,  (string)orank_columns[0]);
   at.set_entry(r, c+1,  (string)orank_columns[1]);
   at.set_entry(r, c+2,  (string)orank_columns[2]);
   at.set_entry(r, c+3,  (string)orank_columns[3]);
   at.set_entry(r, c+4,  (string)orank_columns[4]);
   at.set_entry(r, c+5,  (string)orank_columns[5]);
   at.set_entry(r, c+6,  (string)orank_columns[6]);
   at.set_entry(r, c+7,  (string)orank_columns[7]);
   at.set_entry(r, c+8,  (string)orank_columns[8]);
   at.set_entry(r, c+9,  (string)orank_columns[9]);
   at.set_entry(r, c+10, (string)orank_columns[10]);
   at.set_entry(r, c+11, (string)orank_columns[11]);

   // Write ENS_i for each ensemble member
   for(i=0, col=c+12; i<n_ens; i++) {

     snprintf(tmp_str, sizeof(tmp_str), "%s%i", orank_columns[12], i+1);
      at.set_entry(r, col, (string)tmp_str); // Ensemble member value
      col++;
   }

   at.set_entry(r, c+12+n_ens, (string)orank_columns[13]);
   at.set_entry(r, c+13+n_ens, (string)orank_columns[14]);
   at.set_entry(r, c+14+n_ens, (string)orank_columns[15]);
   at.set_entry(r, c+15+n_ens, (string)orank_columns[16]);
   at.set_entry(r, c+16+n_ens, (string)orank_columns[17]);
   at.set_entry(r, c+17+n_ens, (string)orank_columns[18]);
   at.set_entry(r, c+18+n_ens, (string)orank_columns[19]);
   at.set_entry(r, c+19+n_ens, (string)orank_columns[20]);

   return;
}

////////////////////////////////////////////////////////////////////////

void write_relp_header_row(int hdr_flag, int n_ens, AsciiTable &at,
                           int r, int c) {
   int i, col;
   char tmp_str[max_str_len];

   // Write the header column names if requested
   if(hdr_flag) {
      for(i=0; i<n_header_columns; i++)
         at.set_entry(r, i+c, (string)hdr_columns[i]);

      c += n_header_columns;
   }

   // Write the columns names specific to the RELP line type
   at.set_entry(r, c+0, (string)relp_columns[0]);
   at.set_entry(r, c+1, (string)relp_columns[1]);

   // Write RELP_i for each ensemble member
   for(i=0, col=c+2; i<n_ens; i++) {
     snprintf(tmp_str, sizeof(tmp_str), "%s%i", relp_columns[2], i+1);
      at.set_entry(r, col, (string)tmp_str);
      col++;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_fho_row(StatHdrColumns &shc, const CTSInfo &cts_info,
                   STATOutputType out_type,
                   AsciiTable &stat_at, int &stat_row,
                   AsciiTable &txt_at, int &txt_row) {

   // FHO line type
   shc.set_line_type(stat_fho_str);

   // Thresholds
   shc.set_fcst_thresh(cts_info.fthresh);
   shc.set_obs_thresh(cts_info.othresh);

   // Not Applicable
   shc.set_thresh_logic(SetLogic_None);
   shc.set_alpha(bad_data_double);
   shc.set_cov_thresh(na_str);

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_fho_cols(cts_info, stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   return;
}

////////////////////////////////////////////////////////////////////////

void write_ctc_row(StatHdrColumns &shc, const CTSInfo &cts_info,
                   STATOutputType out_type,
                   AsciiTable &stat_at, int &stat_row,
                   AsciiTable &txt_at, int &txt_row) {

   // CTC line type
   shc.set_line_type(stat_ctc_str);

   // Thresholds
   shc.set_fcst_thresh(cts_info.fthresh);
   shc.set_obs_thresh(cts_info.othresh);

   // Not Applicable
   shc.set_thresh_logic(SetLogic_None);
   shc.set_alpha(bad_data_double);
   shc.set_cov_thresh(na_str);

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_ctc_cols(cts_info, stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   return;
}

////////////////////////////////////////////////////////////////////////

void write_cts_row(StatHdrColumns &shc, const CTSInfo &cts_info,
                   STATOutputType out_type,
                   AsciiTable &stat_at, int &stat_row,
                   AsciiTable &txt_at, int &txt_row) {
   int i;

   // CTS line type
   shc.set_line_type(stat_cts_str);

   // Thresholds
   shc.set_fcst_thresh(cts_info.fthresh);
   shc.set_obs_thresh(cts_info.othresh);

   // Not Applicable
   shc.set_thresh_logic(SetLogic_None);
   shc.set_cov_thresh(na_str);

   // Write a line for each alpha value
   for(i=0; i<cts_info.n_alpha; i++) {

      // Alpha value
      shc.set_alpha(cts_info.alpha[i]);

      // Write the header columns
      write_header_cols(shc, stat_at, stat_row);

      // Write the data columns
      write_cts_cols(cts_info, i, stat_at, stat_row, n_header_columns);

      // If requested, copy row to the text file
      if(out_type == STATOutputType_Both) {
         copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

         // Increment the text row counter
         txt_row++;
      }

      // Increment the STAT row counter
      stat_row++;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_mctc_row(StatHdrColumns &shc, const MCTSInfo &mcts_info,
                    STATOutputType out_type,
                    AsciiTable &stat_at, int &stat_row,
                    AsciiTable &txt_at, int &txt_row) {

   // MCTC line type
   shc.set_line_type(stat_mctc_str);

   // Thresholds
   shc.set_fcst_thresh(mcts_info.fthresh);
   shc.set_obs_thresh(mcts_info.othresh);

   // Not Applicable
   shc.set_thresh_logic(SetLogic_None);
   shc.set_alpha(bad_data_double);
   shc.set_cov_thresh(na_str);

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_mctc_cols(mcts_info, stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   return;
}

////////////////////////////////////////////////////////////////////////

void write_mcts_row(StatHdrColumns &shc, const MCTSInfo &mcts_info,
                    STATOutputType out_type,
                    AsciiTable &stat_at, int &stat_row,
                    AsciiTable &txt_at, int &txt_row) {
   int i;

   // MCTS line type
   shc.set_line_type(stat_mcts_str);

   // Thresholds
   shc.set_fcst_thresh(mcts_info.fthresh);
   shc.set_obs_thresh(mcts_info.othresh);

   // Not Applicable
   shc.set_thresh_logic(SetLogic_None);
   shc.set_cov_thresh(na_str);

   // Write a line for each alpha value
   for(i=0; i<mcts_info.n_alpha; i++) {

      // Alpha value
      shc.set_alpha(mcts_info.alpha[i]);

      // Write the header columns
      write_header_cols(shc, stat_at, stat_row);

      // Write the data columns
      write_mcts_cols(mcts_info, i, stat_at, stat_row, n_header_columns);

      // If requested, copy row to the text file
      if(out_type == STATOutputType_Both) {
         copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

         // Increment the text row counter
         txt_row++;
      }

      // Increment the STAT row counter
      stat_row++;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_cnt_row(StatHdrColumns &shc, const CNTInfo &cnt_info,
                   STATOutputType out_type, int i_bin, int n_bin,
                   AsciiTable &stat_at, int &stat_row,
                   AsciiTable &txt_at, int &txt_row) {
   int i;
   ConcatString mask_name = shc.get_mask();

   // CNT line type
   shc.set_line_type(stat_cnt_str);

   // Thresholds
   shc.set_fcst_thresh(cnt_info.fthresh);
   shc.set_obs_thresh(cnt_info.othresh);
   shc.set_thresh_logic(cnt_info.logic);

   // Not Applicable
   shc.set_cov_thresh(na_str);

   // Update the mask name, if needed.
   ConcatString cs = append_climo_bin(mask_name, i_bin, n_bin);
   shc.set_mask(cs.c_str());

   // Write a line for each alpha value
   for(i=0; i<cnt_info.n_alpha; i++) {

      // Alpha value
      shc.set_alpha(cnt_info.alpha[i]);

      // Write the header columns
      write_header_cols(shc, stat_at, stat_row);

      // Write the data columns
      write_cnt_cols(cnt_info, i, stat_at, stat_row, n_header_columns);

      // If requested, copy row to the text file
      if(out_type == STATOutputType_Both) {
         copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

         // Increment the text row counter
         txt_row++;
      }

      // Increment the STAT row counter
      stat_row++;
   }

   // Reset the mask name
   shc.set_mask(mask_name.c_str());

   return;
}

////////////////////////////////////////////////////////////////////////

void write_sl1l2_row(StatHdrColumns &shc, const SL1L2Info &sl1l2_info,
                     STATOutputType out_type, int i_bin, int n_bin,
                     AsciiTable &stat_at, int &stat_row,
                     AsciiTable &txt_at, int &txt_row) {
   ConcatString mask_name = shc.get_mask();

   // SL1L2 line type
   shc.set_line_type(stat_sl1l2_str);

   // Thresholds
   shc.set_fcst_thresh(sl1l2_info.fthresh);
   shc.set_obs_thresh(sl1l2_info.othresh);
   shc.set_thresh_logic(sl1l2_info.logic);

   // Not Applicable
   shc.set_cov_thresh(na_str);
   shc.set_alpha(bad_data_double);

   // Update the mask name, if needed.
   ConcatString cs = append_climo_bin(mask_name, i_bin, n_bin);
   shc.set_mask(cs.c_str());

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_sl1l2_cols(sl1l2_info, stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   // Reset the mask name
   shc.set_mask(mask_name.c_str());

   return;
}

////////////////////////////////////////////////////////////////////////

void write_sal1l2_row(StatHdrColumns &shc, const SL1L2Info &sl1l2_info,
                      STATOutputType out_type, int i_bin, int n_bin,
                      AsciiTable &stat_at, int &stat_row,
                      AsciiTable &txt_at, int &txt_row) {
   ConcatString mask_name = shc.get_mask();

   // SAL1L2 line type
   shc.set_line_type(stat_sal1l2_str);

   // Thresholds
   shc.set_fcst_thresh(sl1l2_info.fthresh);
   shc.set_obs_thresh(sl1l2_info.othresh);
   shc.set_thresh_logic(sl1l2_info.logic);

   // Not Applicable
   shc.set_cov_thresh(na_str);
   shc.set_alpha(bad_data_double);

   // Update the mask name, if needed.
   ConcatString cs = append_climo_bin(mask_name, i_bin, n_bin);
   shc.set_mask(cs.c_str());

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_sal1l2_cols(sl1l2_info, stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   // Reset the mask name
   shc.set_mask(mask_name.c_str());

   return;
}

////////////////////////////////////////////////////////////////////////

void write_vl1l2_row(StatHdrColumns &shc, const VL1L2Info &vl1l2_info,
                     STATOutputType out_type,
                     AsciiTable &stat_at, int &stat_row,
                     AsciiTable &txt_at, int &txt_row) {

   // VL1L2 line type
   shc.set_line_type(stat_vl1l2_str);

   // Thresholds
   shc.set_fcst_thresh(vl1l2_info.fthresh);
   shc.set_obs_thresh(vl1l2_info.othresh);
   shc.set_thresh_logic(vl1l2_info.logic);

   // Not Applicable
   shc.set_cov_thresh(na_str);
   shc.set_alpha(bad_data_double);

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_vl1l2_cols(vl1l2_info, stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   return;
}

////////////////////////////////////////////////////////////////////////

void write_val1l2_row(StatHdrColumns &shc, const VL1L2Info &vl1l2_info,
                      STATOutputType out_type,
                      AsciiTable &stat_at, int &stat_row,
                      AsciiTable &txt_at, int &txt_row) {

   // VAL1L2 line type
   shc.set_line_type(stat_val1l2_str);

   // Thresholds
   shc.set_fcst_thresh(vl1l2_info.fthresh);
   shc.set_obs_thresh(vl1l2_info.othresh);
   shc.set_thresh_logic(vl1l2_info.logic);

   // Not Applicable
   shc.set_cov_thresh(na_str);
   shc.set_alpha(bad_data_double);

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_val1l2_cols(vl1l2_info, stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   return;
}

////////////////////////////////////////////////////////////////////////


void write_vcnt_row(StatHdrColumns &shc, const VL1L2Info &vcnt_info,
                     STATOutputType out_type,
                     AsciiTable &stat_at, int &stat_row,
                     AsciiTable &txt_at, int &txt_row) {

   // VL1L2 line type
   shc.set_line_type(stat_vcnt_str);

   // Thresholds
   shc.set_fcst_thresh(vcnt_info.fthresh);
   shc.set_obs_thresh(vcnt_info.othresh);
   shc.set_thresh_logic(vcnt_info.logic);

   // Not Applicable
   shc.set_cov_thresh(na_str);
   shc.set_alpha(bad_data_double);

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_vcnt_cols(vcnt_info, stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   return;
}

////////////////////////////////////////////////////////////////////////

void write_pct_row(StatHdrColumns &shc, const PCTInfo &pct_info,
                   STATOutputType out_type, int i_bin, int n_bin,
                   AsciiTable &stat_at, int &stat_row,
                   AsciiTable &txt_at, int &txt_row,
                   bool update_thresh) {
   ConcatString mask_name = shc.get_mask();

   // PCT line type
   shc.set_line_type(stat_pct_str);

   // Set the threshold columns, if requested.
   if(update_thresh) {
      shc.set_fcst_thresh(pct_info.fthresh);
      shc.set_obs_thresh(pct_info.othresh);
      shc.set_thresh_logic(SetLogic_None);
      shc.set_cov_thresh(na_str);
   }

   // Not Applicable
   shc.set_alpha(bad_data_double);

   // Update the mask name, if needed.
   ConcatString cs = append_climo_bin(mask_name, i_bin, n_bin);
   shc.set_mask(cs.c_str());

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_pct_cols(pct_info, stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   // Reset the mask name
   shc.set_mask(mask_name.c_str());

   return;
}

////////////////////////////////////////////////////////////////////////

void write_pstd_row(StatHdrColumns &shc, const PCTInfo &pct_info,
                    STATOutputType out_type, int i_bin, int n_bin,
                    AsciiTable &stat_at, int &stat_row,
                    AsciiTable &txt_at, int &txt_row,
                    bool update_thresh) {
   int i;
   ConcatString mask_name = shc.get_mask();

   // PSTD line type
   shc.set_line_type(stat_pstd_str);

   // Set the threshold columns, if requested.
   if(update_thresh) {
      shc.set_fcst_thresh(pct_info.fthresh);
      shc.set_obs_thresh(pct_info.othresh);
      shc.set_thresh_logic(SetLogic_None);
      shc.set_cov_thresh(na_str);
   }

   // Update the mask name, if needed.
   ConcatString cs = append_climo_bin(mask_name, i_bin, n_bin);
   shc.set_mask(cs.c_str());

   // Write a line for each alpha value
   for(i=0; i<pct_info.n_alpha; i++) {

      // Alpha value
      shc.set_alpha(pct_info.alpha[i]);

      // Write the header columns
      write_header_cols(shc, stat_at, stat_row);

      // Write the data columns
      write_pstd_cols(pct_info, i, stat_at, stat_row, n_header_columns);

      // If requested, copy row to the text file
      if(out_type == STATOutputType_Both) {
         copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

         // Increment the text row counter
         txt_row++;
      }

      // Increment the STAT row counter
      stat_row++;
   }

   // Reset the mask name
   shc.set_mask(mask_name.c_str());

   return;
}

////////////////////////////////////////////////////////////////////////

void write_pjc_row(StatHdrColumns &shc, const PCTInfo &pct_info,
                   STATOutputType out_type, int i_bin, int n_bin,
                   AsciiTable &stat_at, int &stat_row,
                   AsciiTable &txt_at, int &txt_row,
                   bool update_thresh) {
   ConcatString mask_name = shc.get_mask();

   // PJC line type
   shc.set_line_type(stat_pjc_str);

   // Set the threshold columns, if requested.
   if(update_thresh) {
      shc.set_fcst_thresh(pct_info.fthresh);
      shc.set_obs_thresh(pct_info.othresh);
      shc.set_thresh_logic(SetLogic_None);
      shc.set_cov_thresh(na_str);
   }

   // Not Applicable
   shc.set_alpha(bad_data_double);

   // Update the mask name, if needed.
   ConcatString cs = append_climo_bin(mask_name, i_bin, n_bin);
   shc.set_mask(cs.c_str());

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_pjc_cols(pct_info, stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   // Reset the mask name
   shc.set_mask(mask_name.c_str());

   return;
}

////////////////////////////////////////////////////////////////////////

void write_prc_row(StatHdrColumns &shc, const PCTInfo &pct_info,
                   STATOutputType out_type, int i_bin, int n_bin,
                   AsciiTable &stat_at, int &stat_row,
                   AsciiTable &txt_at, int &txt_row,
                   bool update_thresh) {
   ConcatString mask_name = shc.get_mask();

   // PRC line type
   shc.set_line_type(stat_prc_str);

   // Set the threshold columns, if requested.
   if(update_thresh) {
      shc.set_fcst_thresh(pct_info.fthresh);
      shc.set_obs_thresh(pct_info.othresh);
      shc.set_thresh_logic(SetLogic_None);
      shc.set_cov_thresh(na_str);
   }

   // Not Applicable
   shc.set_alpha(bad_data_double);

   // Update the mask name, if needed.
   ConcatString cs = append_climo_bin(mask_name, i_bin, n_bin);
   shc.set_mask(cs.c_str());

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_prc_cols(pct_info, stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   // Reset the mask name
   shc.set_mask(mask_name.c_str());

   return;
}

////////////////////////////////////////////////////////////////////////

void write_eclv_row(StatHdrColumns &shc, const PCTInfo &pct_info,
                    const NumArray &eclv_points,
                    STATOutputType out_type, int i_bin, int n_bin,
                    AsciiTable &stat_at, int &stat_row,
                    AsciiTable &txt_at, int &txt_row) {
   int i;
   ConcatString mask_name = shc.get_mask();

   // ECLV line type
   shc.set_line_type(stat_eclv_str);

   // Set the threshold columns, if requested.
   shc.set_obs_thresh(pct_info.othresh);
   shc.set_thresh_logic(SetLogic_None);
   shc.set_cov_thresh(na_str);

   // Not Applicable
   shc.set_alpha(bad_data_double);

   // Update the mask name, if needed.
   ConcatString cs = append_climo_bin(mask_name, i_bin, n_bin);
   shc.set_mask(cs.c_str());

   // Write ECLV line for each PCT row
   for(i=0; i<pct_info.pct.nrows(); i++) {

      // Update the forecast threshold
      shc.set_fcst_thresh(pct_info.fthresh[i+1]);

      // Write the header columns
      write_header_cols(shc, stat_at, stat_row);

      // Write the data for the 2x2 contingency table for this row
      write_eclv_cols(pct_info.pct.ctc_by_row(i), eclv_points,
                      stat_at, stat_row, n_header_columns);

      // If requested, copy row to the text file
      if(out_type == STATOutputType_Both) {
         copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

         // Increment the text row counter
         txt_row++;
      }

      // Increment the STAT row counter
      stat_row++;

   } // end for i

   // Reset the mask name
   shc.set_mask(mask_name.c_str());

   return;
}

////////////////////////////////////////////////////////////////////////

void write_eclv_row(StatHdrColumns &shc, const CTSInfo &cts_info,
                    const NumArray &eclv_points, STATOutputType out_type,
                    AsciiTable &stat_at, int &stat_row,
                    AsciiTable &txt_at, int &txt_row) {

   // ECLV line type
   shc.set_line_type(stat_eclv_str);

   // Thresholds
   shc.set_fcst_thresh(cts_info.fthresh);
   shc.set_obs_thresh(cts_info.othresh);

   // Not Applicable
   shc.set_thresh_logic(SetLogic_None);
   shc.set_alpha(bad_data_double);
   shc.set_cov_thresh(na_str);

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_eclv_cols(cts_info.cts, eclv_points,
                   stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   return;
}

////////////////////////////////////////////////////////////////////////

void write_nbrctc_row(StatHdrColumns &shc, const NBRCTSInfo &nbrcts_info,
                   STATOutputType out_type,
                   AsciiTable &stat_at, int &stat_row,
                   AsciiTable &txt_at, int &txt_row) {

   // NBRCTC line type
   shc.set_line_type(stat_nbrctc_str);

   // Thresholds
   shc.set_fcst_thresh(nbrcts_info.fthresh);
   shc.set_obs_thresh(nbrcts_info.othresh);

   // Fractional coverage threshold
   shc.set_cov_thresh(nbrcts_info.cthresh);

   // Not Applicable
   shc.set_thresh_logic(SetLogic_None);
   shc.set_alpha(bad_data_double);

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_nbrctc_cols(nbrcts_info, stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   return;
}

////////////////////////////////////////////////////////////////////////

void write_nbrcts_row(StatHdrColumns &shc, const NBRCTSInfo &nbrcts_info,
                      STATOutputType out_type,
                      AsciiTable &stat_at, int &stat_row,
                      AsciiTable &txt_at, int &txt_row) {
   int i;

   // NBRCTS line type
   shc.set_line_type(stat_nbrcts_str);

   // Thresholds
   shc.set_fcst_thresh(nbrcts_info.fthresh);
   shc.set_obs_thresh(nbrcts_info.othresh);

   // Fractional coverage threshold
   shc.set_cov_thresh(nbrcts_info.cthresh);

   // Not Applicable
   shc.set_thresh_logic(SetLogic_None);

   // Write a line for each alpha value
   for(i=0; i<nbrcts_info.cts_info.n_alpha; i++) {

      // Alpha value
      shc.set_alpha(nbrcts_info.cts_info.alpha[i]);

      // Write the header columns
      write_header_cols(shc, stat_at, stat_row);

      // Write the data columns
      write_nbrcts_cols(nbrcts_info, i, stat_at, stat_row,
                        n_header_columns);

      // If requested, copy row to the text file
      if(out_type == STATOutputType_Both) {
         copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

         // Increment the text row counter
         txt_row++;
      }

      // Increment the STAT row counter
      stat_row++;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_nbrcnt_row(StatHdrColumns &shc, const NBRCNTInfo &nbrcnt_info,
                      STATOutputType out_type,
                      AsciiTable &stat_at, int &stat_row,
                      AsciiTable &txt_at, int &txt_row) {
   int i;

   // NBRCNT line type
   shc.set_line_type(stat_nbrcnt_str);

   // Thresholds
   shc.set_fcst_thresh(nbrcnt_info.fthresh);
   shc.set_obs_thresh(nbrcnt_info.othresh);

   // Not applicable
   shc.set_thresh_logic(SetLogic_None);
   shc.set_cov_thresh(na_str);

   // Write a line for each alpha value
   for(i=0; i<nbrcnt_info.n_alpha; i++) {

      // Alpha value
      shc.set_alpha(nbrcnt_info.alpha[i]);

      // Write the header columns
      write_header_cols(shc, stat_at, stat_row);

      // Write the data columns
      write_nbrcnt_cols(nbrcnt_info, i, stat_at, stat_row,
                        n_header_columns);

      // If requested, copy row to the text file
      if(out_type == STATOutputType_Both) {
         copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

         // Increment the text row counter
         txt_row++;
      }

      // Increment the STAT row counter
      stat_row++;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_grad_row(StatHdrColumns &shc, const GRADInfo &grad_info,
                    STATOutputType out_type,
                    AsciiTable &stat_at, int &stat_row,
                    AsciiTable &txt_at, int &txt_row) {

   // GRAD line type
   shc.set_line_type(stat_grad_str);

   // Not applicable
   shc.set_fcst_thresh(na_str);
   shc.set_obs_thresh(na_str);
   shc.set_thresh_logic(SetLogic_None);
   shc.set_cov_thresh(na_str);
   shc.set_alpha(bad_data_double);

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_grad_cols(grad_info, stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   return;
}

////////////////////////////////////////////////////////////////////////

void write_dmap_row(StatHdrColumns &shc, const DMAPInfo &dmap_info,
                    STATOutputType out_type,
                    AsciiTable &stat_at, int &stat_row,
                    AsciiTable &txt_at, int &txt_row) {

   // DMAP line type
   shc.set_line_type(stat_dmap_str);

   // Thresholds
   shc.set_fcst_thresh(dmap_info.fthresh);
   shc.set_obs_thresh(dmap_info.othresh);

   // Not Applicable
   shc.set_thresh_logic(SetLogic_None);
   shc.set_alpha(bad_data_double);
   shc.set_cov_thresh(na_str);

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_dmap_cols(dmap_info, stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   return;
}

////////////////////////////////////////////////////////////////////////

void write_mpr_row(StatHdrColumns &shc, const PairDataPoint *pd_ptr,
                   STATOutputType out_type,
                   AsciiTable &stat_at, int &stat_row,
                   AsciiTable &txt_at, int &txt_row,
                   bool update_thresh) {
   int i;

   // MPR line type
   shc.set_line_type(stat_mpr_str);

   // Set the threshold columns, if requested.
   if(update_thresh) {
      shc.set_fcst_thresh(na_str);
      shc.set_obs_thresh(na_str);
      shc.set_thresh_logic(SetLogic_None);
      shc.set_cov_thresh(na_str);
   }

   // Not Applicable
   shc.set_alpha(bad_data_double);

   // Write a line for each matched pair
   for(i=0; i<pd_ptr->n_obs; i++) {

      // Set the observation valid time
      shc.set_obs_valid_beg(pd_ptr->vld_ta[i]);
      shc.set_obs_valid_end(pd_ptr->vld_ta[i]);

      // Write the header columns
      write_header_cols(shc, stat_at, stat_row);

      // Write the data columns
      write_mpr_cols(pd_ptr, i, stat_at, stat_row, n_header_columns);

      // If requested, copy row to the text file
      if(out_type == STATOutputType_Both) {
         copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

         // Increment the text row counter
         txt_row++;
      }

      // Increment the STAT row counter
      stat_row++;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_isc_row(StatHdrColumns &shc, const ISCInfo &isc_info,
                   STATOutputType out_type,
                   AsciiTable &stat_at, int &stat_row,
                   AsciiTable &txt_at, int &txt_row) {
   int i;

   // ISC line type
   shc.set_line_type(stat_isc_str);

   // Not Applicable
   shc.set_interp_mthd(InterpMthd_None,
                       GridTemplateFactory::GridTemplate_None);
   shc.set_interp_wdth(bad_data_int);
   shc.set_cov_thresh(na_str);
   shc.set_alpha(bad_data_double);

   // Write a line for each scale plus one for the thresholded binary
   // field and one for the father wavelet
   for(i=-1; i<=isc_info.n_scale; i++) {

      // Write the header columns
      write_header_cols(shc, stat_at, stat_row);

      // Write the data columns
      write_isc_cols(isc_info, i, stat_at, stat_row, n_header_columns);

      // If requested, copy row to the text file
      if(out_type == STATOutputType_Both) {
         copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

         // Increment the text row counter
         txt_row++;
      }

      // Increment the STAT row counter
      stat_row++;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_ecnt_row(StatHdrColumns &shc, const ECNTInfo &ecnt_info,
                    STATOutputType out_type,
                    AsciiTable &stat_at, int &stat_row,
                    AsciiTable &txt_at, int &txt_row) {

   // ECNT line type
   shc.set_line_type(stat_ecnt_str);

   // Thresholds
   shc.set_obs_thresh(ecnt_info.othresh);

   // Not Applicable
   shc.set_fcst_thresh(na_str);
   shc.set_thresh_logic(SetLogic_None);
   shc.set_cov_thresh(na_str);
   shc.set_alpha(bad_data_double);

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_ecnt_cols(ecnt_info, stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   return;
}

////////////////////////////////////////////////////////////////////////

void write_rps_row(StatHdrColumns &shc, const RPSInfo &rps_info,
                   STATOutputType out_type,
                   AsciiTable &stat_at, int &stat_row,
                   AsciiTable &txt_at, int &txt_row) {
   ConcatString mask_name = shc.get_mask();

   // RPS line type
   shc.set_line_type(stat_rps_str);

   // Thresholds
   shc.set_fcst_thresh(process_rps_cdp_thresh(rps_info.fthresh));
   shc.set_obs_thresh(rps_info.othresh);

   // Not Applicable
   shc.set_thresh_logic(SetLogic_None);
   shc.set_cov_thresh(na_str);
   shc.set_alpha(bad_data_double);

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_rps_cols(rps_info, stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   return;
}

////////////////////////////////////////////////////////////////////////

void write_rhist_row(StatHdrColumns &shc, const PairDataEnsemble *pd_ptr,
                     STATOutputType out_type,
                     AsciiTable &stat_at, int &stat_row,
                     AsciiTable &txt_at, int &txt_row) {

   // Check for data to write.  Running Ensemble-Stat with skip_const
   // set to true may result in no data.
   if(nint(pd_ptr->rhist_na.sum()) == 0) return;

   // RHIST line type
   shc.set_line_type(stat_rhist_str);

   // Not Applicable
   shc.set_fcst_thresh(na_str);
   shc.set_thresh_logic(SetLogic_None);
   shc.set_cov_thresh(na_str);
   shc.set_alpha(bad_data_double);

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_rhist_cols(pd_ptr, stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   return;
}

////////////////////////////////////////////////////////////////////////

void write_phist_row(StatHdrColumns &shc, const PairDataEnsemble *pd_ptr,
                     STATOutputType out_type,
                     AsciiTable &stat_at, int &stat_row,
                     AsciiTable &txt_at, int &txt_row) {

   // Check for data to write.  Running Ensemble-Stat with skip_const
   // set to true may result in no data.
   if(nint(pd_ptr->phist_na.sum()) == 0) return;

   // PHIST line type
   shc.set_line_type(stat_phist_str);

   // Not Applicable
   shc.set_fcst_thresh(na_str);
   shc.set_thresh_logic(SetLogic_None);
   shc.set_cov_thresh(na_str);
   shc.set_alpha(bad_data_double);

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_phist_cols(pd_ptr, stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   return;
}

////////////////////////////////////////////////////////////////////////

void write_orank_row(StatHdrColumns &shc, const PairDataEnsemble *pd_ptr,
                     STATOutputType out_type,
                     AsciiTable &stat_at, int &stat_row,
                     AsciiTable &txt_at, int &txt_row) {
   int i;

   // Observation Rank line type
   shc.set_line_type(stat_orank_str);

   // Not Applicable
   shc.set_fcst_thresh(na_str);
   shc.set_thresh_logic(SetLogic_None);
   shc.set_cov_thresh(na_str);
   shc.set_alpha(bad_data_double);

   // Write a line for each ensemble pair
   for(i=0; i<pd_ptr->n_obs; i++) {

      // Set the observation valid time
      shc.set_obs_valid_beg(pd_ptr->vld_ta[i]);
      shc.set_obs_valid_end(pd_ptr->vld_ta[i]);

      // Write the header columns
      write_header_cols(shc, stat_at, stat_row);

      // Write the data columns
      write_orank_cols(pd_ptr, i, stat_at, stat_row, n_header_columns);

      // If requested, copy row to the text file
      if(out_type == STATOutputType_Both) {
         copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

         // Increment the text row counter
         txt_row++;
      }

      // Increment the STAT row counter
      stat_row++;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_ssvar_row(StatHdrColumns &shc, const PairDataEnsemble *pd_ptr,
                     double alpha, STATOutputType out_type,
                     AsciiTable &stat_at, int &stat_row,
                     AsciiTable &txt_at, int &txt_row) {
   int i;

   // SSVAR line type
   shc.set_line_type(stat_ssvar_str);

   // Not Applicable
   shc.set_fcst_thresh(na_str);
   shc.set_thresh_logic(SetLogic_None);
   shc.set_cov_thresh(na_str);

   // Alpha value
   shc.set_alpha(alpha);

   // Write a line for each ssvar bin
   for(i=0; i<pd_ptr->ssvar_bins[0].n_bin; i++) {

      // Write the header columns
      write_header_cols(shc, stat_at, stat_row);

      // Write the data columns
      write_ssvar_cols(pd_ptr, i, alpha,
                       stat_at, stat_row, n_header_columns);

      // If requested, copy row to the text file
      if(out_type == STATOutputType_Both) {
         copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

         // Increment the text row counter
         txt_row++;
      }

      // Increment the STAT row counter
      stat_row++;

   } // end for i

   return;
}

////////////////////////////////////////////////////////////////////////

void write_relp_row(StatHdrColumns &shc, const PairDataEnsemble *pd_ptr,
                    STATOutputType out_type,
                    AsciiTable &stat_at, int &stat_row,
                    AsciiTable &txt_at, int &txt_row) {

   // Check for data to write.  Running Ensemble-Stat with skip_const
   // set to true may result in no data.
   if(nint(pd_ptr->relp_na.sum()) == 0) return;

   // RELP line type
   shc.set_line_type(stat_relp_str);

   // Not Applicable
   shc.set_fcst_thresh(na_str);
   shc.set_thresh_logic(SetLogic_None);
   shc.set_cov_thresh(na_str);
   shc.set_alpha(bad_data_double);

   // Write the header columns
   write_header_cols(shc, stat_at, stat_row);

   // Write the data columns
   write_relp_cols(pd_ptr, stat_at, stat_row, n_header_columns);

   // If requested, copy row to the text file
   if(out_type == STATOutputType_Both) {
      copy_ascii_table_row(stat_at, stat_row, txt_at, txt_row);

      // Increment the text row counter
      txt_row++;
   }

   // Increment the STAT row counter
   stat_row++;

   return;
}

////////////////////////////////////////////////////////////////////////

void write_header_cols(const StatHdrColumns &shc,
                       AsciiTable &at, int r) {
   ConcatString cs;

   //
   // Header columns:
   //    VERSION,        MODEL,
   //    DESC,           FCST_LEAD,
   //    FCST_VALID_BEG, FCST_VALID_END,
   //    OBS_LEAD,
   //    OBS_VALID_BEG,  OBS_VALID_END,
   //    FCST_VAR,       FCST_UNITS,
   //    FCST_LEV,       OBS_VAR,
   //    OBS_UNITS,      OBS_LEV,
   //    OBTYPE,         VX_MASK,
   //    INTERP_MTHD,    INTERP_PNTS,
   //    FCST_THRESH,    OBS_THRESH,
   //    COV_THRESH,     ALPHA,
   //    LINE_TYPE
   //
   at.set_entry(r,  0, (string)met_version);          // MET version
   at.set_entry(r,  1, shc.get_model());              // Model name
   at.set_entry(r,  2, shc.get_desc());               // Description
   at.set_entry(r,  3, shc.get_fcst_lead_str());      // Fcst lead time
   at.set_entry(r,  4, shc.get_fcst_valid_beg_str()); // Fcst valid time
   at.set_entry(r,  5, shc.get_fcst_valid_end_str()); // Fcst valid time
   at.set_entry(r,  6, shc.get_obs_lead_str());       // Obs lead time
   at.set_entry(r,  7, shc.get_obs_valid_beg_str());  // Obs valid time
   at.set_entry(r,  8, shc.get_obs_valid_end_str());  // Obs valid time
   at.set_entry(r,  9, shc.get_fcst_var());           // Fcst variable
   at.set_entry(r, 10, shc.get_fcst_units());         // Fcst units
   at.set_entry(r, 11, shc.get_fcst_lev());           // Fcst level
   at.set_entry(r, 12, shc.get_obs_var());            // Obs variable
   at.set_entry(r, 13, shc.get_obs_units());          // Obs units
   at.set_entry(r, 14, shc.get_obs_lev());            // Obs level
   at.set_entry(r, 15, shc.get_obtype());             // Verifying observation type
   at.set_entry(r, 16, shc.get_mask());               // Verification masking region
   at.set_entry(r, 17, shc.get_interp_mthd());        // Interpolation method
   at.set_entry(r, 18, shc.get_interp_pnts_str());    // Interpolation points
   at.set_entry(r, 19, shc.get_fcst_thresh_str());    // Fcst threshold
   at.set_entry(r, 20, shc.get_obs_thresh_str());     // Obs threshold
   at.set_entry(r, 21, shc.get_cov_thresh_str());     // Coverage threshold
   at.set_entry(r, 22, shc.get_alpha());              // Alpha value
   at.set_entry(r, 23, shc.get_line_type());          // Line type

   return;
}

////////////////////////////////////////////////////////////////////////

void write_fho_cols(const CTSInfo &cts_info,
                    AsciiTable &at, int r, int c) {

   //
   // FHO
   // Dump out the FHO line:
   //    TOTAL,       F_RATE,      H_RATE,
   //    O_RATE
   //
   at.set_entry(r, c+0,  // Total Count
      cts_info.cts.n());

   at.set_entry(r, c+1,  // Forecast Rate = FY/N
      cts_info.cts.f_rate());

   at.set_entry(r, c+2,  // Hit Rate = FY_OY/N
      cts_info.cts.h_rate());

   at.set_entry(r, c+3,  // Observation = OY/N
      cts_info.cts.o_rate());

   return;
}

////////////////////////////////////////////////////////////////////////

void write_ctc_cols(const CTSInfo &cts_info,
                    AsciiTable &at, int r, int c) {

   //
   // Contingency Table Counts
   // Dump out the CTC line:
   //    TOTAL,       FY_OY,       FY_ON,
   //    FN_OY,       FN_ON
   //
   at.set_entry(r, c+0,  // Total Count
      cts_info.cts.n());

   at.set_entry(r, c+1,  // FY_OY
      cts_info.cts.fy_oy());

   at.set_entry(r, c+2,  // FY_ON
      cts_info.cts.fy_on());

   at.set_entry(r, c+3,  // FN_OY
      cts_info.cts.fn_oy());

   at.set_entry(r, c+4,  // FN_ON
      cts_info.cts.fn_on());

   return;
}

////////////////////////////////////////////////////////////////////////

void write_cts_cols(const CTSInfo &cts_info, int i,
                    AsciiTable &at, int r, int c) {

   //
   // Contingency Table Stats
   // Dump out the CTS line:
   //    TOTAL,
   //    BASER,       BASER_NCL,   BASER_NCU,   BASER_BCL,   BASER_BCU,
   //    FMEAN,       FMEAN_NCL,   FMEAN_NCU,   FMEAN_BCL,   FMEAN_BCU,
   //    ACC,         ACC_NCL,     ACC_NCU,     ACC_BCL,     ACC_BCU,
   //    FBIAS,       FBIAS_BCL,   FBIAS_BCU,
   //    PODY,        PODY_NCL,    PODY_NCU,    PODY_BCL,    PODY_BCU,
   //    PODN,        PODN_NCL,    PODN_NCU,    PODN_BCL,    PODN_BCU,
   //    POFD,        POFD_NCL,    POFD_NCU,    POFD_BCL,    POFD_BCU,
   //    FAR,         FAR_NCL,     FAR_NCU,     FAR_BCL,     FAR_BCU,
   //    CSI,         CSI_NCL,     CSI_NCU,     CSI_BCL,     CSI_BCU,
   //    GSS,         GSS_BCL,     GSS_BCU,
   //    HK,          HK_NCL,      HK_NCU,      HK_BCL,      HK_BCU,
   //    HSS,         HSS_BCL,     HSS_BCU,
   //    ODDS,        ODDS_NCL,    ODDS_NCU,    ODDS_BCL,    ODDS_BCU,
   //    LOR,         LOR_NCL,     LOR_NCU,     LOR_BCL,     LOR_BCU,
   //    ORSS,        ORSS_NCL,    ORSS_NCU,    ORSS_BCL,    ORSS_BCU,
   //    EDS,         EDS_NCL,     EDS_NCU,     EDS_BCL,     EDS_BCU,
   //    SEDS,        SEDS_NCL,    SEDS_NCU,    SEDS_BCL,    SEDS_BCU,
   //    EDI,         EDI_NCL,     EDI_NCU,     EDI_BCL,     EDI_BCU,
   //    SEDI,        SEDI_NCL,    SEDI_NCU,    SEDI_BCL,    SEDI_BCU,
   //    BAGSS,       BAGSS_BCL,   BAGSS_BCU
   //
   at.set_entry(r, c+0,  // Total count
      cts_info.cts.n());

   at.set_entry(r, c+1,  // Base Rate (oy_tp)
      cts_info.baser.v);

   at.set_entry(r, c+2,  // Base Rate (oy_tp) NCL
      cts_info.baser.v_ncl[i]);

   at.set_entry(r, c+3,  // Base Rate (oy_tp) NCU
      cts_info.baser.v_ncu[i]);

   at.set_entry(r, c+4,  // Base Rate (oy_tp) BCL
      cts_info.baser.v_bcl[i]);

   at.set_entry(r, c+5,  // Base Rate (oy_tp) BCU
      cts_info.baser.v_bcu[i]);

   at.set_entry(r, c+6,  // Forecast Mean (fy_tp)
      cts_info.fmean.v);

   at.set_entry(r, c+7,  // Forecast Mean (fy_tp) NCL
      cts_info.fmean.v_ncl[i]);

   at.set_entry(r, c+8,  // Forecast Mean (fy_tp) NCU
      cts_info.fmean.v_ncu[i]);

   at.set_entry(r, c+9,  // Forecast Mean (fy_tp) BCL
      cts_info.fmean.v_bcl[i]);

   at.set_entry(r, c+10, // Forecast Mean (fy_tp) BCU
      cts_info.fmean.v_bcu[i]);

   at.set_entry(r, c+11, // Accuracy
      cts_info.acc.v);

   at.set_entry(r, c+12, // Accuracy NCL
      cts_info.acc.v_ncl[i]);

   at.set_entry(r, c+13, // Accuracy NCU
      cts_info.acc.v_ncu[i]);

   at.set_entry(r, c+14, // Accuracy BCL
      cts_info.acc.v_bcl[i]);

   at.set_entry(r, c+15, // Accuracy BCU
      cts_info.acc.v_bcu[i]);

   at.set_entry(r, c+16, // Frequency Bias
      cts_info.fbias.v);

   at.set_entry(r, c+17, // Frequency Bias BCL
      cts_info.fbias.v_bcl[i]);

   at.set_entry(r, c+18, // Frequency Bias BCU
      cts_info.fbias.v_bcu[i]);

   at.set_entry(r, c+19, // POD Yes
      cts_info.pody.v);

   at.set_entry(r, c+20, // POD Yes NCL
      cts_info.pody.v_ncl[i]);

   at.set_entry(r, c+21, // POD Yes NCU
      cts_info.pody.v_ncu[i]);

   at.set_entry(r, c+22, // POD Yes BCL
      cts_info.pody.v_bcl[i]);

   at.set_entry(r, c+23, // POD Yes BCU
      cts_info.pody.v_bcu[i]);

   at.set_entry(r, c+24, // POD No
      cts_info.podn.v);

   at.set_entry(r, c+25, // POD No NCL
      cts_info.podn.v_ncl[i]);

   at.set_entry(r, c+26, // POD No NCU
      cts_info.podn.v_ncu[i]);

   at.set_entry(r, c+27, // POD No BCL
      cts_info.podn.v_bcl[i]);

   at.set_entry(r, c+28, // POD No BCU
      cts_info.podn.v_bcu[i]);

   at.set_entry(r, c+29, // POFD
      cts_info.pofd.v);

   at.set_entry(r, c+30, // POFD NCL
      cts_info.pofd.v_ncl[i]);

   at.set_entry(r, c+31, // POFD NCU
      cts_info.pofd.v_ncu[i]);

   at.set_entry(r, c+32, // POFD BCL
      cts_info.pofd.v_bcl[i]);

   at.set_entry(r, c+33, // POFD BCU
      cts_info.pofd.v_bcu[i]);

   at.set_entry(r, c+34, // FAR
      cts_info.far.v);

   at.set_entry(r, c+35, // FAR NCL
      cts_info.far.v_ncl[i]);

   at.set_entry(r, c+36, // FAR NCU
      cts_info.far.v_ncu[i]);

   at.set_entry(r, c+37, // FAR BCL
      cts_info.far.v_bcl[i]);

   at.set_entry(r, c+38, // FAR BCU
      cts_info.far.v_bcu[i]);

   at.set_entry(r, c+39, // CSI (Threat Score)
      cts_info.csi.v);

   at.set_entry(r, c+40, // CSI (Threat Score) NCL
      cts_info.csi.v_ncl[i]);

   at.set_entry(r, c+41, // CSI (Threat Score) NCU
      cts_info.csi.v_ncu[i]);

   at.set_entry(r, c+42, // CSI (Threat Score) BCL
      cts_info.csi.v_bcl[i]);

   at.set_entry(r, c+43, // CSI (Threat Score) BCU
      cts_info.csi.v_bcu[i]);

   at.set_entry(r, c+44, // Gilbert Skill Score (ETS)
      cts_info.gss.v);

   at.set_entry(r, c+45, // Gilbert Skill Score (ETS) BCL
      cts_info.gss.v_bcl[i]);

   at.set_entry(r, c+46, // Gilbert Skill Score (ETS) BCU
      cts_info.gss.v_bcu[i]);

   at.set_entry(r, c+47, // Hanssen-Kuipers Discriminant (TSS)
      cts_info.hk.v);

   at.set_entry(r, c+48, // Hanssen-Kuipers Discriminant (TSS) NCL
      cts_info.hk.v_ncl[i]);

   at.set_entry(r, c+49, // Hanssen-Kuipers Discriminant (TSS) NCU
      cts_info.hk.v_ncu[i]);

   at.set_entry(r, c+50, // Hanssen-Kuipers Discriminant (TSS) BCL
      cts_info.hk.v_bcl[i]);

   at.set_entry(r, c+51, // Hanssen-Kuipers Discriminant (TSS) BCU
      cts_info.hk.v_bcu[i]);

   at.set_entry(r, c+52, // Heidke Skill Score
      cts_info.hss.v);

   at.set_entry(r, c+53, // Heidke Skill Score BCL
      cts_info.hss.v_bcl[i]);

   at.set_entry(r, c+54, // Heidke Skill Score BCU
      cts_info.hss.v_bcu[i]);

   at.set_entry(r, c+55, // Odds Ratio
      cts_info.odds.v);

   at.set_entry(r, c+56, // Odds Ratio NCL
      cts_info.odds.v_ncl[i]);

   at.set_entry(r, c+57, // Odds Ratio NCU
      cts_info.odds.v_ncu[i]);

   at.set_entry(r, c+58, // Odds Ratio BCL
      cts_info.odds.v_bcl[i]);

   at.set_entry(r, c+59, // Odds Ratio BCU
      cts_info.odds.v_bcu[i]);

   at.set_entry(r, c+60, // Log Odds Ratio
      cts_info.lodds.v);

   at.set_entry(r, c+61, // Log Odds Ratio NCL
      cts_info.lodds.v_ncl[i]);

   at.set_entry(r, c+62, // Log Odds Ratio NCU
      cts_info.lodds.v_ncu[i]);

   at.set_entry(r, c+63, // Log Odds Ratio BCL
      cts_info.lodds.v_bcl[i]);

   at.set_entry(r, c+64, // Log Odds Ratio BCU
      cts_info.lodds.v_bcu[i]);

   at.set_entry(r, c+65, // Odds Ratio Skill Score
      cts_info.orss.v);

   at.set_entry(r, c+66, // Odds Ratio Skill Score NCL
      cts_info.orss.v_ncl[i]);

   at.set_entry(r, c+67, // Odds Ratio Skill Score NCU
      cts_info.orss.v_ncu[i]);

   at.set_entry(r, c+68, // Odds Ratio Skill Score BCL
      cts_info.orss.v_bcl[i]);

   at.set_entry(r, c+69, // Odds Ratio Skill Score BCU
      cts_info.orss.v_bcu[i]);

   at.set_entry(r, c+70, // Extreme Dependency Score
      cts_info.eds.v);

   at.set_entry(r, c+71, // Extreme Dependency Score NCL
      cts_info.eds.v_ncl[i]);

   at.set_entry(r, c+72, // Extreme Dependency Score NCU
      cts_info.eds.v_ncu[i]);

   at.set_entry(r, c+73, // Extreme Dependency Score BCL
      cts_info.eds.v_bcl[i]);

   at.set_entry(r, c+74, // Extreme Dependency Score BCU
      cts_info.eds.v_bcu[i]);

   at.set_entry(r, c+75, // Symmetric Extreme Dependency Score
      cts_info.seds.v);

   at.set_entry(r, c+76, // Symmetric Extreme Dependency Score NCL
      cts_info.seds.v_ncl[i]);

   at.set_entry(r, c+77, // Symmetric Extreme Dependency Score NCU
      cts_info.seds.v_ncu[i]);

   at.set_entry(r, c+78, // Symmetric Extreme Dependency Score BCL
      cts_info.seds.v_bcl[i]);

   at.set_entry(r, c+79, // Symmetric Extreme Dependency Score BCU
      cts_info.seds.v_bcu[i]);

   at.set_entry(r, c+80, // Extreme Dependency Index
      cts_info.edi.v);

   at.set_entry(r, c+81, // Extreme Dependency Index NCL
      cts_info.edi.v_ncl[i]);

   at.set_entry(r, c+82, // Extreme Dependency Index NCU
      cts_info.edi.v_ncu[i]);

   at.set_entry(r, c+83, // Extreme Dependency Index BCL
      cts_info.edi.v_bcl[i]);

   at.set_entry(r, c+84, // Extreme Dependency Index BCU
      cts_info.edi.v_bcu[i]);

   at.set_entry(r, c+85, // Symmetric Extreme Dependency Index
      cts_info.sedi.v);

   at.set_entry(r, c+86, // Symmetric Extreme Dependency Index NCL
      cts_info.sedi.v_ncl[i]);

   at.set_entry(r, c+87, // Symmetric Extreme Dependency Index NCU
      cts_info.sedi.v_ncu[i]);

   at.set_entry(r, c+88, // Symmetric Extreme Dependency Index BCL
      cts_info.sedi.v_bcl[i]);

   at.set_entry(r, c+89, // Symmetric Extreme Dependency Index BCU
      cts_info.sedi.v_bcu[i]);

   at.set_entry(r, c+90, // Bias-Corrected Gilbert Skill Score
      cts_info.bagss.v);

   at.set_entry(r, c+91, // Bias-Corrected Gilbert Skill Score BCL
      cts_info.bagss.v_bcl[i]);

   at.set_entry(r, c+92, // Bias-Corrected Gilbert Skill Score BCU
      cts_info.bagss.v_bcu[i]);

   return;
}

////////////////////////////////////////////////////////////////////////

void write_cnt_cols(const CNTInfo &cnt_info, int i,
                    AsciiTable &at, int r, int c) {

   //
   // Continuous Variable Stats
   // Dump out the CNT line:
   //    TOTAL,
   //    FBAR,             FBAR_NCL,             FBAR_NCU,          FBAR_BCL,          FBAR_BCU,
   //    FSTDEV,           FSTDEV_NCL,           FSTDEV_NCU,        FSTDEV_BCL,        FSTDEV_BCU,
   //    OBAR,             OBAR_NCL,             OBAR_NCU,          OBAR_BCL,          OBAR_BCU,
   //    OSTDEV,           OSTDEV_NCL,           OSTDEV_NCU,        OSTDEV_BCL,        OSTDEV_BCU,
   //    PR_CORR,          PR_CORR_NCL,          PR_CORR_NCU,       PR_CORR_BCL,       PR_CORR_BCU,
   //    SP_CORR,          KT_CORR,              RANKS,             FRANK_TIES,        ORANK_TIES,
   //    ME,               ME_NCL,               ME_NCU,            ME_BCL,            ME_BCU,
   //    ESTDEV,           ESTDEV_NCL,           ESTDEV_NCU,        ESTDEV_BCL,        ESTDEV_BCU,
   //    MBIAS,            MBIAS_BCL,            MBIAS_BCU,
   //    MAE,              MAE_BCL,              MAE_BCU,
   //    MSE,              MSE_BCL,              MSE_BCU,
   //    BCMSE,            BCMSE_BCL,            BCMSE_BCU,
   //    RMSE,             RMSE_BCL,             RMSE_BCU,
   //    E10,              E10_BCL,              E10_BCU,
   //    E25,              E25_BCL,              E25_BCU,
   //    E50,              E50_BCL,              E50_BCU,
   //    E75,              E75_BCL,              E75_BCU,
   //    E90,              E90_BCL,              E90_BCU,
   //    EIQR,             EIQR_BCL,             EIQR_BCU,
   //    ANOM_CORR,        ANOM_CORR_NCL,        ANOM_CORR_NCU,     ANOM_CORR_BCL,     ANOM_CORR_BCU,
   //    MAD,              MAD_BCL,              MAD_BCU,
   //    ME2,              ME2_BCL,              ME2_BCU,
   //    MSESS,            MSESS_BCL,            MSESS_BCU,
   //    RMSFA,            RMSFA_BCL,            RMSFA_BCU,
   //    RMSOA,            RMSOA_BCL,            RMSOA_BCU,
   //    ANOM_CORR_UNCNTR, ANOM_CORR_UNCNTR_BCL, ANOM_CORR_UNCNTR_BCU
   //

   at.set_entry(r, c+0,  // Total Number of Grid Points
      cnt_info.n);

   at.set_entry(r, c+1,  // Forecast Mean
      cnt_info.fbar.v);

   at.set_entry(r, c+2,  // Forecast Mean NCL
      cnt_info.fbar.v_ncl[i]);

   at.set_entry(r, c+3,  // Forecast Mean NCU
      cnt_info.fbar.v_ncu[i]);

   at.set_entry(r, c+4,  // Forecast Mean BCL
      cnt_info.fbar.v_bcl[i]);

   at.set_entry(r, c+5,  // Forecast Mean BCU
      cnt_info.fbar.v_bcu[i]);

   at.set_entry(r, c+6,  // Forecast Standard Deviation
      cnt_info.fstdev.v);

   at.set_entry(r, c+7,  // Forecast Standard Deviation NCL
      cnt_info.fstdev.v_ncl[i]);

   at.set_entry(r, c+8,  // Forecast Standard Deviation NCU
      cnt_info.fstdev.v_ncu[i]);

   at.set_entry(r, c+9,  // Forecast Standard Deviation BCL
      cnt_info.fstdev.v_bcl[i]);

   at.set_entry(r, c+10, // Forecast Standard Deviation BCU
      cnt_info.fstdev.v_bcu[i]);

   at.set_entry(r, c+11, // Observation Mean
      cnt_info.obar.v);

   at.set_entry(r, c+12, // Observation Mean NCL
      cnt_info.obar.v_ncl[i]);

   at.set_entry(r, c+13, // Observation Mean NCU
      cnt_info.obar.v_ncu[i]);

   at.set_entry(r, c+14, // Observation Mean BCL
      cnt_info.obar.v_bcl[i]);

   at.set_entry(r, c+15, // Observation Mean BCU
      cnt_info.obar.v_bcu[i]);

   at.set_entry(r, c+16, // Observation Standard Deviation
      cnt_info.ostdev.v);

   at.set_entry(r, c+17, // Observation Standard Deviation NCL
      cnt_info.ostdev.v_ncl[i]);

   at.set_entry(r, c+18, // Observation Standard Deviation NCU
      cnt_info.ostdev.v_ncu[i]);

   at.set_entry(r, c+19, // Observation Standard Deviation BCL
      cnt_info.ostdev.v_bcl[i]);

   at.set_entry(r, c+20, // Observation Standard Deviation BCU
      cnt_info.ostdev.v_bcu[i]);

   at.set_entry(r, c+21, // Pearson's Correlation Coefficient
      cnt_info.pr_corr.v);

   at.set_entry(r, c+22, // Pearson's Correlation Coefficient NCL
      cnt_info.pr_corr.v_ncl[i]);

   at.set_entry(r, c+23, // Pearson's Correlation Coefficient NCU
      cnt_info.pr_corr.v_ncu[i]);

   at.set_entry(r, c+24, // Pearson's Correlation Coefficient BCL
      cnt_info.pr_corr.v_bcl[i]);

   at.set_entry(r, c+25, // Pearson's Correlation Coefficient BCU
      cnt_info.pr_corr.v_bcu[i]);

   at.set_entry(r, c+26, // Spearman's Rank Correlation Coefficient
      cnt_info.sp_corr.v);

   at.set_entry(r, c+27, // Kendall Tau Rank Correlation Coefficient
      cnt_info.kt_corr.v);

   at.set_entry(r, c+28, // Number of ranks compared
      cnt_info.n_ranks);

   at.set_entry(r, c+29, // Number of tied forecast ranks
      cnt_info.frank_ties);

   at.set_entry(r, c+30, // Number of tied observation ranks
      cnt_info.orank_ties);

   at.set_entry(r, c+31, // Mean Error
      cnt_info.me.v);

   at.set_entry(r, c+32, // Mean Error NCL
      cnt_info.me.v_ncl[i]);

   at.set_entry(r, c+33, // Mean Error NCU
      cnt_info.me.v_ncu[i]);

   at.set_entry(r, c+34, // Mean Error BCL
      cnt_info.me.v_bcl[i]);

   at.set_entry(r, c+35, // Mean Error BCU
      cnt_info.me.v_bcu[i]);

   at.set_entry(r, c+36, // Error Standard Deviation
      cnt_info.estdev.v);

   at.set_entry(r, c+37, // Error Standard Deviation NCL
      cnt_info.estdev.v_ncl[i]);

   at.set_entry(r, c+38, // Error Standard Deviation NCU
      cnt_info.estdev.v_ncu[i]);

   at.set_entry(r, c+39, // Error Standard Deviation BCL
      cnt_info.estdev.v_bcl[i]);

   at.set_entry(r, c+40, // Error Standard Deviation BCU
      cnt_info.estdev.v_bcu[i]);

   at.set_entry(r, c+41, // Multiplicative Bias
      cnt_info.mbias.v);

   at.set_entry(r, c+42, // Multiplicative Bias BCL
      cnt_info.mbias.v_bcl[i]);

   at.set_entry(r, c+43, // Multiplicative Bias BCU
      cnt_info.mbias.v_bcu[i]);

   at.set_entry(r, c+44, // Mean Absolute Error
      cnt_info.mae.v);

   at.set_entry(r, c+45, // Mean Absolute Error BCL
      cnt_info.mae.v_bcl[i]);

   at.set_entry(r, c+46, // Mean Absolute Error BCU
      cnt_info.mae.v_bcu[i]);

   at.set_entry(r, c+47, // Mean Squared Error
      cnt_info.mse.v);

   at.set_entry(r, c+48, // Mean Squared Error BCL
      cnt_info.mse.v_bcl[i]);

   at.set_entry(r, c+49, // Mean Squared Error BCU
      cnt_info.mse.v_bcu[i]);

   at.set_entry(r, c+50, // Bias-Corrected Mean Squared Error
      cnt_info.bcmse.v);

   at.set_entry(r, c+51, // Bias-Corrected Mean Squared Error BCL
      cnt_info.bcmse.v_bcl[i]);

   at.set_entry(r, c+52, // Bias-Corrected Mean Squared Error BCU
      cnt_info.bcmse.v_bcu[i]);

   at.set_entry(r, c+53, // Root Mean Squared Error
      cnt_info.rmse.v);

   at.set_entry(r, c+54, // Root Mean Squared Error BCL
      cnt_info.rmse.v_bcl[i]);

   at.set_entry(r, c+55, // Root Mean Squared Error BCU
      cnt_info.rmse.v_bcu[i]);

   at.set_entry(r, c+56, // 10th Percentile of the Error
      cnt_info.e10.v);

   at.set_entry(r, c+57, // 10th Percentile of the Error BCL
      cnt_info.e10.v_bcl[i]);

   at.set_entry(r, c+58, // 10th Percentile of the Error BCU
      cnt_info.e10.v_bcu[i]);

   at.set_entry(r, c+59, // 25th Percentile of the Error
      cnt_info.e25.v);

   at.set_entry(r, c+60, // 25th Percentile of the Error BCL
      cnt_info.e25.v_bcl[i]);

   at.set_entry(r, c+61, // 25th Percentile of the Error BCU
      cnt_info.e25.v_bcu[i]);

   at.set_entry(r, c+62, // 50th Percentile of the Error
      cnt_info.e50.v);

   at.set_entry(r, c+63, // 50th Percentile of the Error BCL
      cnt_info.e50.v_bcl[i]);

   at.set_entry(r, c+64, // 50th Percentile of the Error BCU
      cnt_info.e50.v_bcu[i]);

   at.set_entry(r, c+65, // 75th Percentile of the Error
      cnt_info.e75.v);

   at.set_entry(r, c+66, // 75th Percentile of the Error BCL
      cnt_info.e75.v_bcl[i]);

   at.set_entry(r, c+67, // 75th Percentile of the Error BCU
      cnt_info.e75.v_bcu[i]);

   at.set_entry(r, c+68, // 90th Percentile of the Error
      cnt_info.e90.v);

   at.set_entry(r, c+69, // 90th Percentile of the Error BCL
      cnt_info.e90.v_bcl[i]);

   at.set_entry(r, c+70, // 90th Percentile of the Error BCU
      cnt_info.e90.v_bcu[i]);

   at.set_entry(r, c+71, // Interquartile Range of the Error
      cnt_info.eiqr.v);

   at.set_entry(r, c+72, // Interquartile Range of the Error BCL
      cnt_info.eiqr.v_bcl[i]);

   at.set_entry(r, c+73, // Interquartile Range of the Error BCU
      cnt_info.eiqr.v_bcu[i]);

   at.set_entry(r, c+74, // Median Absolute Deviation
      cnt_info.mad.v);

   at.set_entry(r, c+75, // Median Absolute Deviation BCL
      cnt_info.mad.v_bcl[i]);

   at.set_entry(r, c+76, // Median Absolute Deviation BCU
      cnt_info.mad.v_bcu[i]);

   at.set_entry(r, c+77, // Anomaly Correlation
      cnt_info.anom_corr.v);

   at.set_entry(r, c+78, // Anomaly Correlation NCL
      cnt_info.anom_corr.v_ncl[i]);

   at.set_entry(r, c+79, // Anomaly Correlation NCU
      cnt_info.anom_corr.v_ncu[i]);

   at.set_entry(r, c+80, // Anomaly Correlation BCL
      cnt_info.anom_corr.v_bcl[i]);

   at.set_entry(r, c+81, // Anomaly Correlation BCU
      cnt_info.anom_corr.v_bcu[i]);

   at.set_entry(r, c+82, // Mean Error Squared
      cnt_info.me2.v);

   at.set_entry(r, c+83, // Mean Error Squared BCL
      cnt_info.me2.v_bcl[i]);

   at.set_entry(r, c+84, // Mean Error Squared BCU
      cnt_info.me2.v_bcu[i]);

   at.set_entry(r, c+85, // Mean Squared Error Skill Score
      cnt_info.msess.v);

   at.set_entry(r, c+86, // Mean Squared Error Skill Score BCL
      cnt_info.msess.v_bcl[i]);

   at.set_entry(r, c+87, // Mean Squared Error Skill Score BCU
      cnt_info.msess.v_bcu[i]);

   at.set_entry(r, c+88, // Root Mean Squared Forecast Anomaly
      cnt_info.rmsfa.v);

   at.set_entry(r, c+89, // Root Mean Squared Forecast Anomaly BCL
      cnt_info.rmsfa.v_bcl[i]);

   at.set_entry(r, c+90, // Root Mean Squared Forecast Anomaly BCU
      cnt_info.rmsfa.v_bcu[i]);

   at.set_entry(r, c+91, // Root Mean Squared Observation Anomaly
      cnt_info.rmsoa.v);

   at.set_entry(r, c+92, // Root Mean Squared Observation Anomaly BCL
      cnt_info.rmsoa.v_bcl[i]);

   at.set_entry(r, c+93, // Root Mean Squared Observation Anomaly BCU
      cnt_info.rmsoa.v_bcu[i]);

   at.set_entry(r, c+94, // Anomaly Correlation Uncentered
      cnt_info.anom_corr_uncntr.v);

   at.set_entry(r, c+95, // Anomaly Correlation Uncentered BCL
      cnt_info.anom_corr_uncntr.v_bcl[i]);

   at.set_entry(r, c+96, // Anomaly Correlation Uncentered BCU
      cnt_info.anom_corr_uncntr.v_bcu[i]);

   return;
}

////////////////////////////////////////////////////////////////////////

void write_mctc_cols(const MCTSInfo &mcts_info,
                     AsciiTable &at, int r, int c) {
   int i, j, col;

   //
   // Multi-Category Contingency Table Counts
   // Dump out the MCTC line:
   //    TOTAL,       N_CAT,     Fi_Oj,     EC_VALUE
   //
   at.set_entry(r, c+0,  // Total Count
      mcts_info.cts.total());

   at.set_entry(r, c+1,  // Number of categories 
      mcts_info.cts.nrows());

   //
   // Loop through the contingency table rows and columns
   //
   for(i=0, col=c+2; i<mcts_info.cts.nrows(); i++) {
      for(j=0; j<mcts_info.cts.ncols(); j++) {

         at.set_entry(r, col,      // Fi_Oj table counts
            mcts_info.cts.entry(i, j));
         col++;
      }
   }

   at.set_entry(r, col,  // Expected Correct value 
      mcts_info.cts.ec_value());
   col++;

   return;
}

////////////////////////////////////////////////////////////////////////

void write_mcts_cols(const MCTSInfo &mcts_info, int i,
                     AsciiTable &at, int r, int c) {

   //
   // Multi-Category Contingency Table Stats
   // Dump out the MCTS line:
   //    TOTAL,       N_CAT,
   //    ACC,         ACC_NCL,     ACC_NCU,     ACC_BCL,     ACC_BCU,
   //    HK,          HK_BCL,      HK_BCU,
   //    HSS,         HSS_BCL,     HSS_BCU,
   //    GER,         GER_BCL,     GER_BCU,
   //    HSS_EC,      HSS_EC_BCL,  HSS_EC_BCU,
   //    EC_VALUE
   //
   at.set_entry(r, c+0,  // Total count
      mcts_info.cts.total());

   at.set_entry(r, c+1,  // Number of categories
      mcts_info.cts.nrows());

   at.set_entry(r, c+2, // Accuracy
      mcts_info.acc.v);

   at.set_entry(r, c+3, // Accuracy NCL
      mcts_info.acc.v_ncl[i]);

   at.set_entry(r, c+4, // Accuracy NCU
      mcts_info.acc.v_ncu[i]);

   at.set_entry(r, c+5, // Accuracy BCL
      mcts_info.acc.v_bcl[i]);

   at.set_entry(r, c+6, // Accuracy BCU
      mcts_info.acc.v_bcu[i]);

   at.set_entry(r, c+7, // Hanssen-Kuipers Discriminant (TSS)
      mcts_info.hk.v);

   at.set_entry(r, c+8, // Hanssen-Kuipers Discriminant (TSS) BCL
      mcts_info.hk.v_bcl[i]);

   at.set_entry(r, c+9, // Hanssen-Kuipers Discriminant (TSS) BCU
      mcts_info.hk.v_bcu[i]);

   at.set_entry(r, c+10, // Heidke Skill Score
      mcts_info.hss.v);

   at.set_entry(r, c+11, // Heidke Skill Score BCL
      mcts_info.hss.v_bcl[i]);

   at.set_entry(r, c+12, // Heidke Skill Score BCU
      mcts_info.hss.v_bcu[i]);

   at.set_entry(r, c+13, // Gerrity Score
      mcts_info.ger.v);

   at.set_entry(r, c+14, // Gerrity Score BCL
      mcts_info.ger.v_bcl[i]);

   at.set_entry(r, c+15, // Gerrity Score BCU
      mcts_info.ger.v_bcu[i]);

   at.set_entry(r, c+16, // Heidke Skill Score with Expected Correct
      mcts_info.hss_ec.v);

   at.set_entry(r, c+17, // Heidke Skill Score EC BCL
      mcts_info.hss_ec.v_bcl[i]);

   at.set_entry(r, c+18, // Heidke Skill Score EC BCU
      mcts_info.hss_ec.v_bcu[i]);

   at.set_entry(r, c+19, // Expected Correct value
      mcts_info.cts.ec_value());

   return;
}

////////////////////////////////////////////////////////////////////////

void write_sl1l2_cols(const SL1L2Info &sl1l2_info,
                      AsciiTable &at, int r, int c) {

   //
   // Scalar L1L2 Line Type (SL1L2)
   // Dump out the SL1L2 line:
   //    TOTAL,       FBAR,        OBAR,
   //    FOBAR,       FFBAR,       OOBAR
   //
   at.set_entry(r, c+0,  // Total Count
      sl1l2_info.scount);

   at.set_entry(r, c+1,  // FBAR
      sl1l2_info.fbar);

   at.set_entry(r, c+2,  // OBAR
      sl1l2_info.obar);

   at.set_entry(r, c+3,  // FOBAR
      sl1l2_info.fobar);

   at.set_entry(r, c+4,  // FFBAR
      sl1l2_info.ffbar);

   at.set_entry(r, c+5,  // OOBAR
      sl1l2_info.oobar);

   at.set_entry(r, c+6,  // MAE
      sl1l2_info.mae);

   return;
}

////////////////////////////////////////////////////////////////////////

void write_sal1l2_cols(const SL1L2Info &sl1l2_info,
                       AsciiTable &at, int r, int c) {

   //
   // Scalar Anomaly L1L2 Line Type (SAL1L2)
   // Dump out the SAL1L2 line:
   //    TOTAL,       FABAR,       OABAR,
   //    FOABAR,      FFABAR,      OOABAR,
   //    MAE
   //
   at.set_entry(r, c+0,  // Total Anomaly Count
      sl1l2_info.sacount);

   at.set_entry(r, c+1,  // FABAR
      sl1l2_info.fabar);

   at.set_entry(r, c+2,  // OABAR
      sl1l2_info.oabar);

   at.set_entry(r, c+3,  // FOABAR
      sl1l2_info.foabar);

   at.set_entry(r, c+4,  // FFABAR
      sl1l2_info.ffabar);

   at.set_entry(r, c+5,  // OOABAR
      sl1l2_info.ooabar);

   at.set_entry(r, c+6,  // MAE
      sl1l2_info.mae);

   return;
}

////////////////////////////////////////////////////////////////////////

void write_vl1l2_cols(const VL1L2Info &vl1l2_info, AsciiTable &at, int r, int c)

{

   //
   // Vector L1L2 Line Type (VL1L2)
   // Dump out the VL1L2 line:
   //    TOTAL,       UFBAR,       VFBAR,
   //    UOBAR,       VOBAR,       UVFOBAR,
   //    UVFFBAR,     UVOOBAR      F_SPEED_BAR,
   //    O_SPEED_BAR,
   //

   at.set_entry(r, c+0,  // Total Count
      vl1l2_info.vcount);

   at.set_entry(r, c+1,  // UFBAR
      vl1l2_info.uf_bar);

   at.set_entry(r, c+2,  // VFBAR
      vl1l2_info.vf_bar);

   at.set_entry(r, c+3,  // UOBAR
      vl1l2_info.uo_bar);

   at.set_entry(r, c+4,  // VOBAR
      vl1l2_info.vo_bar);

   at.set_entry(r, c+5,  // UVFOBAR
      vl1l2_info.uvfo_bar);

   at.set_entry(r, c+6,  // UVFFBAR
      vl1l2_info.uvff_bar);

   at.set_entry(r, c+7,  // UVOOBAR
      vl1l2_info.uvoo_bar);

   at.set_entry(r, c+8,  // F_SPEED_BAR
      vl1l2_info.f_speed_bar);

   at.set_entry(r, c+9,  // O_SPEED_BAR
      vl1l2_info.o_speed_bar);

   return;
}

////////////////////////////////////////////////////////////////////////

void write_val1l2_cols(const VL1L2Info &vl1l2_info, AsciiTable &at, int r, int c)

{

   //
   // Vector Anomaly L1L2 Line Type (VAL1L2)
   // Dump out the VAL1L2 line:
   //    TOTAL,       UFABAR,      VFABAR,
   //    UOABAR,      VOABAR,      UVFOABAR,
   //    UVFFABAR,    UVOOABAR
   //

   at.set_entry(r, c+0,  // Total Anomaly Count
      vl1l2_info.vacount);

   at.set_entry(r, c+1,  // UFABAR
      vl1l2_info.ufa_bar);

   at.set_entry(r, c+2,  // VFABAR
      vl1l2_info.vfa_bar);

   at.set_entry(r, c+3,  // UOABAR
      vl1l2_info.uoa_bar);

   at.set_entry(r, c+4,  // VOABAR
      vl1l2_info.voa_bar);

   at.set_entry(r, c+5,  // UVFOABAR
      vl1l2_info.uvfoa_bar);

   at.set_entry(r, c+6,  // UVFFABAR
      vl1l2_info.uvffa_bar);

   at.set_entry(r, c+7,  // UVOOABAR
      vl1l2_info.uvooa_bar);
   return;
}

////////////////////////////////////////////////////////////////////////


void write_vcnt_cols(const VL1L2Info &vcnt_info, AsciiTable &at, int r, int c)

{

   //
   // VCNT Line Type
   //

     // TOTAL,
     //
     // FBAR,         FBAR_BCL,        FBAR_BCU,
     // OBAR,         OBAR_BCL,        OBAR_BCU,
     // FS_RMS,       FS_RMS_BCL,      FS_RMS_BCU,
     // OS_RMS,       OS_RMS_BCL,      OS_RMS_BCU,
     // MSVE,         MSVE_BCL,        MSVE_BCU,
     // RMSVE,        RMSVE_BCL,       RMSVE_BCU,
     // FSTDEV,       FSTDEV_BCL,      FSTDEV_BCU,
     // OSTDEV,       OSTDEV_BCL,      OSTDEV_BCU,
     // FDIR,         FDIR_BCL,        FDIR_BCU,
     // ODIR,         ODIR_BCL,        ODIR_BCU,
     // FBAR_SPEED,   FBAR_SPEED_BCL,  FBAR_SPEED_BCU,
     // OBAR_SPEED,   OBAR_SPEED_BCL,  OBAR_SPEED_BCU,
     // VDIFF_SPEED,  VDIFF_SPEED_BCL, VDIFF_SPEED_BCU,
     // VDIFF_DIR,    VDIFF_DIR_BCL,   VDIFF_DIR_BCU,
     // SPD_ERR,      SPD_ERR_BCL,     SPD_ERR_BCU,
     // SPD_ABSERR,   SPD_ABSERR_BCL,  SPD_ABSERR_BCU,
     // DIR_ERR,      DIR_ERR_BCL,     DIR_ERR_BCU,
     // DIR_ABSERR,   DIR_ABSERR_BCL,  DIR_ABSERR_BCU,


   at.set_entry(r, c++, vcnt_info.vcount);         // TOTAL

   at.set_entry(r, c++, vcnt_info.FBAR);           // FBAR
   at.set_entry(r, c++, (string)na_str);                   // FBAR_BCL
   at.set_entry(r, c++, (string)na_str);                   // FBAR_BCU

   at.set_entry(r, c++, vcnt_info.OBAR);           // OBAR
   at.set_entry(r, c++, (string)na_str);                   // OBAR_BCL
   at.set_entry(r, c++, (string)na_str);                   // OBAR_BCU

   at.set_entry(r, c++, vcnt_info.FS_RMS);         // FS_RMS
   at.set_entry(r, c++, (string)na_str);                   // FS_RMS_BCL
   at.set_entry(r, c++, (string)na_str);                   // FS_RMS_BCU

   at.set_entry(r, c++, vcnt_info.OS_RMS);         // OS_RMS
   at.set_entry(r, c++, (string)na_str);                   // OS_RMS_BCL
   at.set_entry(r, c++, (string)na_str);                   // OS_RMS_BCU

   at.set_entry(r, c++, vcnt_info.MSVE);           // MSVE
   at.set_entry(r, c++, (string)na_str);                   // MSVE_BCL
   at.set_entry(r, c++, (string)na_str);                   // MSVE_BCU

   at.set_entry(r, c++, vcnt_info.RMSVE);          // RMSVE
   at.set_entry(r, c++, (string)na_str);                   // RMSVE_BCL
   at.set_entry(r, c++, (string)na_str);                   // RMSVE_BCU

   at.set_entry(r, c++, vcnt_info.FSTDEV);         // FSTDEV
   at.set_entry(r, c++, (string)na_str);                   // FSTDEV_BCL
   at.set_entry(r, c++, (string)na_str);                   // FSTDEV_BCU

   at.set_entry(r, c++, vcnt_info.OSTDEV);         // OSTDEV
   at.set_entry(r, c++, (string)na_str);                   // OSTDEV_BCL
   at.set_entry(r, c++, (string)na_str);                   // OSTDEV_BCU

   at.set_entry(r, c++, vcnt_info.FDIR);           // FDIR
   at.set_entry(r, c++, (string)na_str);                   // FDIR_BCL
   at.set_entry(r, c++, (string)na_str);                   // FDIR_BCU

   at.set_entry(r, c++, vcnt_info.ODIR);           // ODIR
   at.set_entry(r, c++, (string)na_str);                   // ODIR_BCL
   at.set_entry(r, c++, (string)na_str);                   // ODIR_BCU

   at.set_entry(r, c++, vcnt_info.FBAR_SPEED);     // FBAR_SPEED
   at.set_entry(r, c++, (string)na_str);                   // FBAR_SPEED_BCL
   at.set_entry(r, c++, (string)na_str);                   // FBAR_SPEED_BCU

   at.set_entry(r, c++, vcnt_info.OBAR_SPEED);     // OBAR_SPEED
   at.set_entry(r, c++, (string)na_str);                   // OBAR_SPEED_BCL
   at.set_entry(r, c++, (string)na_str);                   // OBAR_SPEED_BCU

   at.set_entry(r, c++, vcnt_info.VDIFF_SPEED);    // VDIFF_SPEED
   at.set_entry(r, c++, (string)na_str);                   // VDIFF_SPEED_BCL
   at.set_entry(r, c++, (string)na_str);                   // VDIFF_SPEED_BCU

   at.set_entry(r, c++, vcnt_info.VDIFF_DIR);      // VDIFF_DIR
   at.set_entry(r, c++, (string)na_str);                   // VDIFF_DIR_BCL
   at.set_entry(r, c++, (string)na_str);                   // VDIFF_DIR_BCU

   at.set_entry(r, c++, vcnt_info.SPEED_ERR);      // SPEED_ERR
   at.set_entry(r, c++, (string)na_str);                   // SPEED_ERR_BCL
   at.set_entry(r, c++, (string)na_str);                   // SPEED_ERR_BCU

   at.set_entry(r, c++, vcnt_info.SPEED_ABSERR);   // SPEED_ABSERR
   at.set_entry(r, c++, (string)na_str);                   // SPEED_ABSERR_BCL
   at.set_entry(r, c++, (string)na_str);                   // SPEED_ABSERR_BCU

   at.set_entry(r, c++, vcnt_info.DIR_ERR);        // DIR_ERR
   at.set_entry(r, c++, (string)na_str);                   // DIR_ERR_BCL
   at.set_entry(r, c++, (string)na_str);                   // DIR_ERR_BCU

   at.set_entry(r, c++, vcnt_info.DIR_ABSERR);     // DIR_ABSERR
   at.set_entry(r, c++, (string)na_str);                   // DIR_ABSERR_BCL
   at.set_entry(r, c++, (string)na_str);                   // DIR_ABSERR_BCU

   //
   //
   //

return;

}


////////////////////////////////////////////////////////////////////////

void write_pct_cols(const PCTInfo &pct_info,
                    AsciiTable &at, int r, int c) {
   int i, col;

   //
   // Nx2 Contingency Table Counts for Probability Forecast
   // Dump out the PCT line:
   //    TOTAL,       N_THRESH,
   //   [THRESH,      OY,          ON,] (for each row of Nx2 Table)
   //    THRESH                         (last threshold)
   //
   at.set_entry(r, c+0,    // Total Count
      pct_info.pct.n());

   at.set_entry(r, c+1,    // N_THRESH
      pct_info.pct.nrows() + 1);

   //
   // Write THRESH_i, OY_i, ON_i for each row of the Nx2 table
   //
   for(i=0, col=c+2; i<pct_info.pct.nrows(); i++) {

      at.set_entry(r, col, // THRESH
         pct_info.pct.threshold(i));
      col++;

      at.set_entry(r, col, // Event Count (OY)
         pct_info.pct.event_count_by_row(i));
      col++;

      at.set_entry(r, col, // Non-Event Count (ON)
         pct_info.pct.nonevent_count_by_row(i));
      col++;
   }

   //
   // Write out the last threshold
   //
   at.set_entry(r, col,    // THRESH
      pct_info.pct.threshold(pct_info.pct.nrows()));

   return;
}

////////////////////////////////////////////////////////////////////////

void write_pstd_cols(const PCTInfo &pct_info, int alpha_i,
                     AsciiTable &at, int r, int c) {
   int i, col;

   //
   // Nx2 Contingency Table Statistics for Probability Forecast
   // Dump out the PSTD line:
   //    TOTAL,       N_THRESH,    BASER,
   //    BASER_NCL,   BASER_NCU,   RELIABILTY,
   //    RESOLUTION,  UNCERTAINTY, ROC_AUC,
   //    BRIER,       BRIER_NCL,   BRIER_NCU,
   //    BRIERCL,     BRIERCL_NCL, BRIERCL_NCU,
   //    BSS,         BSS_SMPL,    [THRESH] (for each threshold)
   //
   at.set_entry(r, c+0,  // Total count
      pct_info.total);

   at.set_entry(r, c+1,  // N_THRESH
      pct_info.pct.nrows() + 1);

   at.set_entry(r, c+2,  // BASER
      pct_info.baser.v);

   at.set_entry(r, c+3,  // BASER_NCL
      pct_info.baser.v_ncl[alpha_i]);

   at.set_entry(r, c+4,  // BASER_NCU
      pct_info.baser.v_ncu[alpha_i]);

   at.set_entry(r, c+5,  // RELIABILITY
      pct_info.reliability);

   at.set_entry(r, c+6,  // RESOLUTION
      pct_info.resolution);

   at.set_entry(r, c+7,  // UNCERTAINTY
      pct_info.uncertainty);

   at.set_entry(r, c+8,  // ROC_AUC
      pct_info.roc_auc);

   at.set_entry(r, c+9,  // BRIER
      pct_info.brier.v);

   at.set_entry(r, c+10, // BRIER_NCL
      pct_info.brier.v_ncl[alpha_i]);

   at.set_entry(r, c+11, // BRIER_NCU
      pct_info.brier.v_ncu[alpha_i]);

   at.set_entry(r, c+12,  // BRIERCL
      pct_info.briercl.v);

   at.set_entry(r, c+13, // BRIERCL_NCL
      pct_info.briercl.v_ncl[alpha_i]);

   at.set_entry(r, c+14, // BRIERCL_NCU
      pct_info.briercl.v_ncu[alpha_i]);

   at.set_entry(r, c+15, // Brier Skill Score
      pct_info.bss);

   at.set_entry(r, c+16, // Brier Skill Score using sample climo
      pct_info.bss_smpl);

   //
   // Write THRESH_i for each probability threshold
   //
   for(i=0, col=c+17; i<=pct_info.pct.nrows(); i++) {

      at.set_entry(r, col, // THRESH
         pct_info.pct.threshold(i));
      col++;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_pjc_cols(const PCTInfo &pct_info,
                    AsciiTable &at, int r, int c) {
   int i, col, n;

   //
   // Nx2 Contingency Table Joint/Continuous Probability
   // Dump out the PJC line:
   //    TOTAL,       N_THRESH,
   //   [THRESH,
   //    OY_TP,       ON_TP,      CALIBRATION,
   //    REFINEMENT,  LIKELIHOOD, BASER] (for each row of Nx2 Table)
   //    THRESH                          (last threshold)
   //
   at.set_entry(r, c+0,    // Total Count
      pct_info.pct.n());

   at.set_entry(r, c+1,    // N_THRESH
      pct_info.pct.nrows() + 1);

   //
   // Write THRESH, OY, ON for each row of the Nx2 table
   //
   n = pct_info.pct.n();
   for(i=0, col=c+2; i<pct_info.pct.nrows(); i++) {

      at.set_entry(r, col, // THRESH
         pct_info.pct.threshold(i));
      col++;

      at.set_entry(r, col, // OY_TP
         pct_info.pct.event_count_by_row(i)/(double) n);
      col++;

      at.set_entry(r, col, // ON_TP
         pct_info.pct.nonevent_count_by_row(i)/(double) n);
      col++;

      at.set_entry(r, col, // CALIBRATION
         pct_info.pct.row_calibration(i));
      col++;

      at.set_entry(r, col, // REFINEMENT
         pct_info.pct.row_refinement(i));
      col++;

      at.set_entry(r, col, // LIKELIHOOD
         pct_info.pct.row_event_likelihood(i));
      col++;

      at.set_entry(r, col, // BASER
         pct_info.pct.row_obar(i));
      col++;
   }

   //
   // Write out the last threshold
   //
   at.set_entry(r, col,    // THRESH
      pct_info.pct.threshold(pct_info.pct.nrows()));

   return;
}

////////////////////////////////////////////////////////////////////////

void write_prc_cols(const PCTInfo &pct_info,
                    AsciiTable &at, int r, int c) {
   int i, col;
   TTContingencyTable ct;

   //
   // Nx2 Contingency Table Points for the ROC Curve
   // Dump out the PRC line:
   //    TOTAL,       N_THRESH,
   //   [THRESH,      PODY,        POFD,] (for each row of Nx2 Table)
   //    THRESH                           (last threshold)
   //
   at.set_entry(r, c+0,    // Total Count
      pct_info.pct.n());

   at.set_entry(r, c+1,    // N_THRESH
      pct_info.pct.nrows() + 1);

   //
   // Write THRESH, PODY, POFD for each row of the Nx2 table
   //
   for(i=0, col=c+2; i<pct_info.pct.nrows(); i++) {

      //
      // Get the 2x2 contingency table for this row
      //
      ct = pct_info.pct.ctc_by_row(i);

      at.set_entry(r, col, // THRESH
         pct_info.pct.threshold(i));
      col++;

      at.set_entry(r, col, // PODY
         ct.pod_yes());
      col++;

      at.set_entry(r, col, // POFD
         ct.pofd());
      col++;
   }

   //
   // Write out the last threshold
   //
   at.set_entry(r, col,    // THRESH
      pct_info.pct.threshold(pct_info.pct.nrows()));

   return;
}

////////////////////////////////////////////////////////////////////////

void write_eclv_cols(const TTContingencyTable &ct,
                     const NumArray &eclv_points,
                     AsciiTable &at, int r, int c) {
   int i, col;

   //
   // Economic Cost/Loss Value
   // Dump out the ECLV line:
   //    TOTAL,   BASER,   BASER_VALUE,
   //    N_PNT,   [CL_],   [VALUE_] (for each point)
   //
   at.set_entry(r, c+0,  // Total Number of pairs
      ct.n());

   at.set_entry(r, c+1,  // Base Rate
      ct.baser());

   at.set_entry(r, c+2,  // Cost/Loss value for Base Rate
      ct.cost_loss(ct.baser()));

   at.set_entry(r, c+3,  // Number of points
      eclv_points.n_elements());

   //
   // Write CL_i and VALUE_i count for each bin
   //
   for(i=0, col=c+4; i<eclv_points.n_elements(); i++) {

      at.set_entry(r, col, // CL_i
         eclv_points[i]);
      col++;

      at.set_entry(r, col, // VALUE_i
         ct.cost_loss(eclv_points[i]));
      col++;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_nbrctc_cols(const NBRCTSInfo &nbrcts_info,
                       AsciiTable &at, int r, int c) {

   //
   // Neighborhood Method Contingency Table Counts Line Type (NBRCTC)
   //
   write_ctc_cols(nbrcts_info.cts_info, at, r, c);

   return;
}

////////////////////////////////////////////////////////////////////////

void write_nbrcts_cols(const NBRCTSInfo &nbrcts_info, int i,
                       AsciiTable &at, int r, int c) {

   //
   // Neighborhood Method Contingency Table Statistics Line Type (NBRCTS)
   //
   write_cts_cols(nbrcts_info.cts_info, i, at, r, c);

   return;
}

////////////////////////////////////////////////////////////////////////

void write_nbrcnt_cols(const NBRCNTInfo &nbrcnt_info, int i,
                       AsciiTable &at, int r, int c) {

   //
   // Neighborhood Method Continuous Statistics (NBRCNT)
   // Dump out the NBRCNT line:
   //    TOTAL,
   //    FBS,         FBS_BCL,     FBS_BCU,
   //    FSS,         FSS_BCL,     FSS_BCU,
   //    AFSS,        AFSS_BCL,    AFSS_BCU,
   //    UFSS,        UFSS_BCL,    UFSS_BCU,
   //    F_RATE,      F_RATE_BCL,  F_RATE_BCU,
   //    O_RATE,      O_RATE_BCL,  O_RATE_BCU
   //
   at.set_entry(r, c+0,  // Total Count
      nbrcnt_info.sl1l2_info.scount);

   at.set_entry(r, c+1,  // Fractions Brier Score
      nbrcnt_info.fbs.v);

   at.set_entry(r, c+2,  // Fractions Brier Score BCL
      nbrcnt_info.fbs.v_bcl[i]);

   at.set_entry(r, c+3,  // Fractions Brier Score BCU
      nbrcnt_info.fbs.v_bcu[i]);

   at.set_entry(r, c+4,  // Fractions Skill Score
      nbrcnt_info.fss.v);

   at.set_entry(r, c+5,  // Fractions Skill Score BCL
      nbrcnt_info.fss.v_bcl[i]);

   at.set_entry(r, c+6,  // Fractions Skill Score BCU
      nbrcnt_info.fss.v_bcu[i]);

   at.set_entry(r, c+7,  // Asymptotic Fractions Skill Score
      nbrcnt_info.afss.v);

   at.set_entry(r, c+8,  // Asymptotic Fractions Skill Score BCL
      nbrcnt_info.afss.v_bcl[i]);

   at.set_entry(r, c+9,  // Asymptotic Fractions Skill Score BCU
      nbrcnt_info.afss.v_bcu[i]);

   at.set_entry(r, c+10, // Uniform Fractions Skill Score
      nbrcnt_info.ufss.v);

   at.set_entry(r, c+11, // Uniform Fractions Skill Score BCL
      nbrcnt_info.ufss.v_bcl[i]);

   at.set_entry(r, c+12, // Uniform Fractions Skill Score BCU
      nbrcnt_info.ufss.v_bcu[i]);

   at.set_entry(r, c+13, // Forecast Rate
      nbrcnt_info.f_rate.v);

   at.set_entry(r, c+14, // Forecast Rate BCL
      nbrcnt_info.f_rate.v_bcl[i]);

   at.set_entry(r, c+15, // Forecast Rate BCU
      nbrcnt_info.f_rate.v_bcu[i]);

   at.set_entry(r, c+16, // Observation Rate
      nbrcnt_info.o_rate.v);

   at.set_entry(r, c+17, // Observation Rate BCL
      nbrcnt_info.o_rate.v_bcl[i]);

   at.set_entry(r, c+18, // Observation Rate BCU
      nbrcnt_info.o_rate.v_bcu[i]);

   return;
}

////////////////////////////////////////////////////////////////////////

void write_grad_cols(const GRADInfo &grad_info,
                     AsciiTable &at, int r, int c) {

   //
   // Gradient Line Type (GRAD)
   //    TOTAL,
   //    FGBAR,       OGBAR,       MGBAR,
   //    EGBAR,       S1,          S1_OG,
   //    FGOG_RATIO,  DX,          DY
   //
   at.set_entry(r, c+0,  // Total Count
      grad_info.total);

   at.set_entry(r, c+1,  // FGBAR
      grad_info.fgbar);

   at.set_entry(r, c+2,  // OGBAR
      grad_info.ogbar);

   at.set_entry(r, c+3,  // MGBAR
      grad_info.mgbar);

   at.set_entry(r, c+4,  // EGBAR
      grad_info.egbar);

   at.set_entry(r, c+5,  // S1
      grad_info.s1());

   at.set_entry(r, c+6,  // S1_OG
      grad_info.s1_og());

   at.set_entry(r, c+7,  // FGOG_RATIO
      grad_info.fgog_ratio());

   at.set_entry(r, c+8,  // DX
      grad_info.dx);

   at.set_entry(r, c+9,  // DY
      grad_info.dy);

   return;
}

////////////////////////////////////////////////////////////////////////

void write_dmap_cols(const DMAPInfo &dmap_info,
                     AsciiTable &at, int r, int c) {

   //
   // Distance Map Line Type (DMAP)
   //    TOTAL,       FY,          OY,
   //    FBIAS,       BADDELEY,    HAUSDORFF,
   //    MED_FO,      MED_OF,      MED_MIN,      MED_MAX,      MED_MEAN,
   //    FOM_FO,      FOM_OF,      FOM_MIN,      FOM_MAX,      FOM_MEAN,
   //    ZHU_FO,      ZHU_OF,      ZHU_MIN,      ZHU_MAX,      ZHU_MEAN,
   //    G,           GBETA,       BETA_VALUE
   //
   at.set_entry(r, c+0,  // TOTAL
      dmap_info.total);

   at.set_entry(r, c+1,  // FY
      dmap_info.fy);

   at.set_entry(r, c+2,  // OY
      dmap_info.oy);

   at.set_entry(r, c+3,  // FBIAS
      dmap_info.fbias());

   at.set_entry(r, c+4,  // BADDELEY
      dmap_info.baddeley);

   at.set_entry(r, c+5,  // HAUSDORFF
      dmap_info.hausdorff);

   at.set_entry(r, c+6,  // MED_FO
      dmap_info.med_fo);

   at.set_entry(r, c+7,  // MED_OF
      dmap_info.med_of);

   at.set_entry(r, c+8,  // MED_MIN
      dmap_info.med_min);

   at.set_entry(r, c+9,  // MED_MAX
      dmap_info.med_max);

   at.set_entry(r, c+10, // MED_MEAN
      dmap_info.med_mean);

   at.set_entry(r, c+11, // FOM_FO
      dmap_info.fom_fo);

   at.set_entry(r, c+12, // FOM_OF
      dmap_info.fom_of);

   at.set_entry(r, c+13, // FOM_MIN
      dmap_info.fom_min);

   at.set_entry(r, c+14, // FOM_MAX
      dmap_info.fom_max);

   at.set_entry(r, c+15, // FOM_MEAN
      dmap_info.fom_mean);

   at.set_entry(r, c+16, // ZHU_FO
      dmap_info.zhu_fo);

   at.set_entry(r, c+17, // ZHU_OF
      dmap_info.zhu_of);

   at.set_entry(r, c+18, // ZHU_MIN
      dmap_info.zhu_min);

   at.set_entry(r, c+19, // ZHU_MAX
      dmap_info.zhu_max);

   at.set_entry(r, c+20, // ZHU_MEAN
      dmap_info.zhu_mean);

   at.set_entry(r, c+21, // G
      dmap_info.g);

   at.set_entry(r, c+22, // GBETA
      dmap_info.gbeta);

   at.set_entry(r, c+23, // BETA_VALUE 
      dmap_info.get_beta_value());

   return;
}

////////////////////////////////////////////////////////////////////////

void write_mpr_cols(const PairDataPoint *pd_ptr, int i,
                    AsciiTable &at, int r, int c) {

   //
   // Matched Pairs (MPR)
   // Dump out the MPR line:
   //    TOTAL,       INDEX,       OBS_SID,
   //    OBS_LAT,     OBS_LON,     OBS_LVL,
   //    OBS_ELV,     FCST,        OBS,
   //    OBS_QC,      CLIMO_MEAN,  CLIMO_STDEV,
   //    CLIMO_CDF
   //
   at.set_entry(r, c+0,  // Total Number of Pairs
      pd_ptr->n_obs);

   at.set_entry(r, c+1,  // Index of Current Pair
      i+1);

   at.set_entry(r, c+2,  // Station ID
      (string)pd_ptr->sid_sa[i]);

   at.set_entry(r, c+3,  // Latitude
      pd_ptr->lat_na[i]);

   at.set_entry(r, c+4,  // Longitude
      pd_ptr->lon_na[i]);

   at.set_entry(r, c+5,  // Level
      pd_ptr->lvl_na[i]);

   at.set_entry(r, c+6,  // Elevation
      pd_ptr->elv_na[i]);

   at.set_entry(r, c+7,  // Forecast Value
      pd_ptr->f_na[i]);

   at.set_entry(r, c+8,  // Observation Value
      pd_ptr->o_na[i]);

   at.set_entry(r, c+9,  // Observation Quality Control
      (string)pd_ptr->o_qc_sa[i]);

   at.set_entry(r, c+10, // Climatological Mean Value
      pd_ptr->cmn_na[i]);

   at.set_entry(r, c+11, // Climatological Standard Deviation Value
      pd_ptr->csd_na[i]);

   at.set_entry(r, c+12, // Climatological CDF Value
      pd_ptr->cdf_na[i]);

   return;
}

////////////////////////////////////////////////////////////////////////

void write_isc_cols(const ISCInfo &isc_info, int i,
                    AsciiTable &at, int r, int c) {

   //
   // Intensity Scale Statistics (ISC)
   // Dump out the ISC line:
   //    TOTAL,
   //    NSCALE,      ISCALE,      MSE,
   //    ISC,         FENERGY2,    OENERGY2,
   //    BASER,       FBIAS
   //
   at.set_entry(r, c+0,     // Total Number of Points
      isc_info.total);

   at.set_entry(r, c+1,     // Tile dimension
      isc_info.tile_dim);

   at.set_entry(r, c+2,     // Tile x-lower left point
      isc_info.tile_xll);

   at.set_entry(r, c+3,     // Tile y-lower left point
      isc_info.tile_yll);

   at.set_entry(r, c+4,     // Total Number of Scales
      isc_info.n_scale);

   at.set_entry(r, c+5,     // Index of Current Scale
      i+1);

   if(i < 0) {
      at.set_entry(r, c+6,  // Mean Squared Error
         isc_info.mse);

      at.set_entry(r, c+7,  // Intensity Scale Score
         isc_info.isc);

      at.set_entry(r, c+8,  // Forecast Energy Squared
         isc_info.fen);

      at.set_entry(r, c+9,  // Observation Energy Squared
         isc_info.oen);
   }
   else {
      at.set_entry(r, c+6,  // Mean Squared Error
         isc_info.mse_scale[i]);

      at.set_entry(r, c+7,  // Intensity Scale Score
         isc_info.isc_scale[i]);

      at.set_entry(r, c+8,  // Forecast Energy Squared
         isc_info.fen_scale[i]);

      at.set_entry(r, c+9,  // Observation Energy Squared
         isc_info.oen_scale[i]);
   }

   at.set_entry(r, c+10,    // Base Rate
      isc_info.baser);

   at.set_entry(r, c+11,    // Forecast Bias
      isc_info.fbias);

   return;
}

////////////////////////////////////////////////////////////////////////

void write_ecnt_cols(const ECNTInfo &ecnt_info,
                     AsciiTable &at, int r, int c) {

   //
   // Ensemble Continuous Statistics
   // Dump out the ECNT line:
   //    TOTAL,            N_ENS,        CRPS,
   //    CRPSS,            IGN,          ME,
   //    RMSE,             SPREAD,       ME_OERR,
   //    RMSE_OERR,        SPREAD_OERR,  SPREAD_PLUS_OERR,
   //    CRPSCL,           CRPS_EMP,     CRPSCL_EMP,
   //    CRPSS_EMP
   //
   at.set_entry(r, c+0,  // Total Number of Pairs
      ecnt_info.n_pair);

   at.set_entry(r, c+1,  // Number of ensemble members
      ecnt_info.n_ens);

   at.set_entry(r, c+2,  // Continuous Ranked Probability Score
      ecnt_info.crps_gaus);

   at.set_entry(r, c+3,  // Continuous Ranked Probability Skill Score
      ecnt_info.crpss_gaus);

   at.set_entry(r, c+4,  // Ignorance Score
      ecnt_info.ign);

   at.set_entry(r, c+5,  // ME for unperturbed ensemble mean
      ecnt_info.me);

   at.set_entry(r, c+6,  // RMSE for unperturbed ensemble mean
      ecnt_info.rmse);

   at.set_entry(r, c+7,  // Mean of unperturbed ensemble spread
      ecnt_info.spread);

   at.set_entry(r, c+8,  // ME for mean of perturbed members
      ecnt_info.me_oerr);

   at.set_entry(r, c+9,  // RMSE for mean of perturbed members
      ecnt_info.rmse_oerr);

   at.set_entry(r, c+10,  // Mean of perturbed ensemble spread
      ecnt_info.spread_oerr);

   at.set_entry(r, c+11,  // Mean of unperturbed spread plus observation error
      ecnt_info.spread_plus_oerr);

   at.set_entry(r, c+12,  // Gaussian climatological CRPS
      ecnt_info.crpscl_gaus);

   at.set_entry(r, c+13,  // Empirical ensemble CRPS
      ecnt_info.crps_emp);

   at.set_entry(r, c+14,  // Empirical climatological CRPS
      ecnt_info.crpscl_emp);

   at.set_entry(r, c+15,  // Empirical CRPSS
      ecnt_info.crpss_emp);

   return;
}

////////////////////////////////////////////////////////////////////////

void write_rps_cols(const RPSInfo &rps_info,
                    AsciiTable &at, int r, int c) {

   //
   // Ensemble Continuous Statistics
   // Dump out the RPS line:
   //    TOTAL,        N_PROB,
   //    RPS_REL,      RPS_RES,    RPS_UNC,
   //    RPS,          RPSS,       RPSS_SMPL,
   //    RPS_COMP
   //
   at.set_entry(r, c+0,  // Total Number of Pairs
      rps_info.n_pair);

   at.set_entry(r, c+1,  // Number of probability bins
      rps_info.n_prob);

   at.set_entry(r, c+2,  // RPS Reliability
      rps_info.rps_rel);

   at.set_entry(r, c+3,  // RPS Resolution
      rps_info.rps_res);

   at.set_entry(r, c+4,  // RPS Uncertainty
      rps_info.rps_unc);

   at.set_entry(r, c+5,  // Ranked Probability Score
      rps_info.rps);

   at.set_entry(r, c+6,  // Ranked Probability Skill Score
      rps_info.rpss);

   at.set_entry(r, c+7,  // Ranked Probability Score using sample climo
      rps_info.rpss_smpl);

   at.set_entry(r, c+8,  // Complement of the Ranked Probability Score
      rps_info.rps_comp());

   return;
}

////////////////////////////////////////////////////////////////////////

void write_rhist_cols(const PairDataEnsemble *pd_ptr,
                      AsciiTable &at, int r, int c) {
   int i, col;

   //
   // Ensemble Ranked Histogram
   // Dump out the RHIST line:
   //    TOTAL,   N_RANK,
   //    [RANK_] (for each bin)
   //
   at.set_entry(r, c+0,  // Total Number of Ranked Observations
      nint(pd_ptr->rhist_na.sum()));

   at.set_entry(r, c+1,  // Total Number of Ranks
      pd_ptr->rhist_na.n_elements());

   //
   // Write RANK_i count for each bin
   //
   for(i=0, col=c+2; i<pd_ptr->rhist_na.n_elements(); i++) {

      at.set_entry(r, col, // RANK_i
         nint(pd_ptr->rhist_na[i]));
      col++;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_phist_cols(const PairDataEnsemble *pd_ptr,
                      AsciiTable &at, int r, int c) {
   int i, col;

   //
   // Probability Integral Transform Histogram
   // Dump out the PHIST line:
   //    TOTAL, BIN_SIZE, N_BIN, [BIN_] (for each bin)
   //
   at.set_entry(r, c+0,  // Total Number of PIT values
      nint(pd_ptr->phist_na.sum()));

   at.set_entry(r, c+1,  // Bin size
      pd_ptr->phist_bin_size);

   at.set_entry(r, c+2,  // Number of bins
      pd_ptr->phist_na.n_elements());

   //
   // Write BIN_i count for each bin
   //
   for(i=0, col=c+3; i<pd_ptr->phist_na.n_elements(); i++) {

      at.set_entry(r, col, // BIN_i
         nint(pd_ptr->phist_na[i]));
      col++;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void write_orank_cols(const PairDataEnsemble *pd_ptr, int i,
                      AsciiTable &at, int r, int c) {
   int j, col;

   //
   // Ensemble Observation Rank Matched Pairs
   // Dump out the ORANK line:
   //    TOTAL,            INDEX,         OBS_SID,
   //    OBS_LAT,          OBS_LON,       OBS_LVL,
   //    OBS_ELV,          OBS,           PIT,
   //    RANK,             N_ENS_VLD,     N_ENS,
   //    [ENS_] (for each ensemble member)
   //    OBS_QC,           ENS_MEAN,      CLIMO_MEAN,
   //    SPREAD,           ENS_MEAN_OERR, SPREAD_OERR,
   //    SPREAD_PLUS_OERR, CLIMO_STDEV
   //
   at.set_entry(r, c+0,  // Total Number of Pairs
      pd_ptr->n_obs);    // Use n_obs instead of n_pair to include missing data

   at.set_entry(r, c+1,  // Index of Current Pair
      i+1);

   at.set_entry(r, c+2,  // Station ID
      (string)pd_ptr->sid_sa[i]);

   at.set_entry(r, c+3,  // Latitude
      pd_ptr->lat_na[i]);

   at.set_entry(r, c+4,  // Longitude
      pd_ptr->lon_na[i]);

   at.set_entry(r, c+5,  // Level
      pd_ptr->lvl_na[i]);

   at.set_entry(r, c+6,  // Elevation
      pd_ptr->elv_na[i]);

   at.set_entry(r, c+7,  // Observation Value
      pd_ptr->o_na[i]);

   at.set_entry(r, c+8,  // Probability Integral Transform
      pd_ptr->pit_na[i]);

   at.set_entry(r, c+9,  // Observation Rank
      nint(pd_ptr->r_na[i]));

   at.set_entry(r, c+10, // Number of valid ensembles
      nint(pd_ptr->v_na[i]));

   at.set_entry(r, c+11, // Number of ensembles
      pd_ptr->n_ens);

   //
   // Write ENS_j for each ensemble member
   //
   for(j=0, col=c+12; j<pd_ptr->n_ens; j++) {

      at.set_entry(r, col, // ENS_j
         pd_ptr->e_na[j][i]);
      col++;
   }

   // Observation Quality Control
   at.set_entry(r, c+12+pd_ptr->n_ens,
      (string)pd_ptr->o_qc_sa[i]);

   // Unperturbed ensemble mean values
   at.set_entry(r, c+13+pd_ptr->n_ens,
      pd_ptr->mn_na[i]);

   // Climatology mean values
   at.set_entry(r, c+14+pd_ptr->n_ens,
      pd_ptr->cmn_na[i]);

   // Unperturbed ensemble spread values
   at.set_entry(r, c+15+pd_ptr->n_ens,
      square_root(pd_ptr->var_na[i]));

   // Perturbed ensemble mean values
   at.set_entry(r, c+16+pd_ptr->n_ens,
      pd_ptr->mn_oerr_na[i]);

   // Perturbed ensemble spread values
   at.set_entry(r, c+17+pd_ptr->n_ens,
      square_root(pd_ptr->var_oerr_na[i]));

   // Unperturbed ensemble spread values plus observation error
   at.set_entry(r, c+18+pd_ptr->n_ens,
      square_root(pd_ptr->var_plus_oerr_na[i]));

   // Climatology standard deviation values
   at.set_entry(r, c+19+pd_ptr->n_ens,
      pd_ptr->csd_na[i]);

   return;
}

////////////////////////////////////////////////////////////////////////

void write_ssvar_cols(const PairDataEnsemble *pd_ptr, int i,
                      double alpha, AsciiTable &at, int r, int c) {
   CNTInfo cnt_info;

   //
   // Allocate space for confidence intervals and derive continuous
   // statistics
   //
   cnt_info.allocate_n_alpha(1);
   cnt_info.alpha[0] = alpha;
   compute_cntinfo(pd_ptr->ssvar_bins[i].sl1l2_info, 0, cnt_info);

   //
   // Ensemble spread/skill variance bins
   // Dump out the SSVAR line:
   //    TOTAL,       N_BIN,       BIN_i,
   //    BIN_N,       VAR_MIN,     VAR_MAX,
   //    VAR_MEAN,    FBAR,        OBAR,
   //    FOBAR,       FFBAR,       OOBAR,
   //    FBAR_NCL,    FBAR_NCU,
   //    FSTDEV,      FSTDEV_NCL,  FSTDEV_NCU,
   //    OBAR_NCL,    OBAR_NCU,
   //    OSTDEV,      OSTDEV_NCL,  OSTDEV_NCU,
   //    PR_CORR,     PR_CORR_NCL, PR_CORR_NCU,
   //    ME,          ME_NCL,      ME_NCU,
   //    ESTDEV,      ESTDEV_NCL,  ESTDEV_NCU,
   //    MBIAS,       MSE,         BCMSE,
   //    RMSE
   //

   at.set_entry(r, c+0,  // Total Number of Pairs
      pd_ptr->n_pair);

   at.set_entry(r, c+1,  // Total Number of Bins
      pd_ptr->ssvar_bins[i].n_bin);

   at.set_entry(r, c+2,  // Index of current bin
      pd_ptr->ssvar_bins[i].bin_i);

   at.set_entry(r, c+3,  // Number of points in bin i
      pd_ptr->ssvar_bins[i].bin_n);

   at.set_entry(r, c+4,  // Lower variance value for bin
      pd_ptr->ssvar_bins[i].var_min);

   at.set_entry(r, c+5,  // Upper variance value for bin
      pd_ptr->ssvar_bins[i].var_max);

   at.set_entry(r, c+6,  // Mean variance value for bin
      pd_ptr->ssvar_bins[i].var_mean);

   at.set_entry(r, c+7,  // Forecast Mean
      pd_ptr->ssvar_bins[i].sl1l2_info.fbar);

   at.set_entry(r, c+8,  // Observation Mean
      pd_ptr->ssvar_bins[i].sl1l2_info.obar);

   at.set_entry(r, c+9,  // Forecast times Observation Mean
      pd_ptr->ssvar_bins[i].sl1l2_info.fobar);

   at.set_entry(r, c+10, // Forecast Squared Mean
      pd_ptr->ssvar_bins[i].sl1l2_info.ffbar);

   at.set_entry(r, c+11, // Observation Squared Mean
      pd_ptr->ssvar_bins[i].sl1l2_info.oobar);

   at.set_entry(r, c+12, // Forecast Mean NCL
      cnt_info.fbar.v_ncl[0]);

   at.set_entry(r, c+13, // Forecast Mean NCU
      cnt_info.fbar.v_ncu[0]);

   at.set_entry(r, c+14, // Forecast Standard Deviation
      cnt_info.fstdev.v);

   at.set_entry(r, c+15, // Forecast Standard Deviation NCL
      cnt_info.fstdev.v_ncl[0]);

   at.set_entry(r, c+16, // Forecast Standard Deviation NCU
      cnt_info.fstdev.v_ncu[0]);

   at.set_entry(r, c+17, // Observation Mean NCL
      cnt_info.obar.v_ncl[0]);

   at.set_entry(r, c+18, // Observation Mean NCU
      cnt_info.obar.v_ncu[0]);

   at.set_entry(r, c+19, // Observation Standard Deviation
      cnt_info.ostdev.v);

   at.set_entry(r, c+20, // Observation Standard Deviation NCL
      cnt_info.ostdev.v_ncl[0]);

   at.set_entry(r, c+21, // Observation Standard Deviation NCU
      cnt_info.ostdev.v_ncu[0]);

   at.set_entry(r, c+22, // Pearson's Correlation Coefficient
      cnt_info.pr_corr.v);

   at.set_entry(r, c+23, // Pearson's Correlation Coefficient NCL
      cnt_info.pr_corr.v_ncl[0]);

   at.set_entry(r, c+24, // Pearson's Correlation Coefficient NCU
      cnt_info.pr_corr.v_ncu[0]);

   at.set_entry(r, c+25, // Mean Error
      cnt_info.me.v);

   at.set_entry(r, c+26, // Mean Error NCL
      cnt_info.me.v_ncl[0]);

   at.set_entry(r, c+27, // Mean Error NCU
      cnt_info.me.v_ncu[0]);

   at.set_entry(r, c+28, // Error Standard Deviation
      cnt_info.estdev.v);

   at.set_entry(r, c+29, // Error Standard Deviation NCL
      cnt_info.estdev.v_ncl[0]);

   at.set_entry(r, c+30, // Error Standard Deviation NCU
      cnt_info.estdev.v_ncu[0]);

   at.set_entry(r, c+31, // Multiplicative Bias
      cnt_info.mbias.v);

   at.set_entry(r, c+32, // Mean Squared Error
      cnt_info.mse.v);

   at.set_entry(r, c+33, // Bias-Corrected Mean Squared Error
      cnt_info.bcmse.v);

   at.set_entry(r, c+34, // Root Mean Squared Error
      cnt_info.rmse.v);

   return;
}

////////////////////////////////////////////////////////////////////////

void write_relp_cols(const PairDataEnsemble *pd_ptr,
                     AsciiTable &at, int r, int c) {
   int i, col;

   //
   // Relative Position
   // Dump out the RELP line:
   //    TOTAL,   N_ENS, [RELP_] (for each ensemble member)
   //
   at.set_entry(r, c+0,  // Total Number of pairs
      nint(pd_ptr->relp_na.sum()));

   at.set_entry(r, c+1,  // Ensemble size
      pd_ptr->relp_na.n_elements());

   //
   // Write RELP_i count for each bin
   //
   for(i=0, col=c+2; i<pd_ptr->relp_na.n_elements(); i++) {

      at.set_entry(r, col, // RELP_i
         pd_ptr->relp_na[i]);
      col++;
   }

   return;
}

////////////////////////////////////////////////////////////////////////

void justify_stat_cols(AsciiTable &at) {

   justify_met_at(at, n_header_columns);

   return;
}

////////////////////////////////////////////////////////////////////////
