.plot_h2 <- function(x) {
    # x is likely the output of compare_h2 (data.frame), NOT fa_asreml
    # But plot() dispatch expects fa_asreml or the user calls plot(h2_results) directly?
    # Wait, plot.fa_asreml is for fa_asreml objects. compare_h2 returns a dataframe.
    # User needs to call plot() on the *dataframe* if they want to plot h2 results?
    # Or we store h2 results inside fa_asreml? No, compare_h2 is separate.

    # Correction: The previous design in roadmap said:
    # "Function: plot(fa_obj, type = 'h2')"
    # But fa_obj doesn't store H2. H2 is computed by compare_h2.
    # So either we compute H2 on the fly inside plot.fa_asreml (requires model, can't do),
    # OR we define plot.data.frame (risky), OR we define a new class for compare_h2 results.

    # Decision: define class "h2_comparison" for the output of compare_h2.
    # But for now, since I'm editing plot.fa_asreml, let's assume the USER passes the dataframe
    # directly to a new plot method, or I add a new S3 method.

    # However, I already added 'type="h2"' to plot.fa_asreml. This implies I should extract something from fa_asreml.
    # But Vg is in fa_asreml. PEV is not.
    # So 'type="h2"' inside plot.fa_asreml is impossible without the model.

    # Re-plan: The user should plot the result of compare_h2.
    # I will revert the change to plot.fa_asreml and instead add a plot method for 'h2_comparison' class.
    # I need to update compare_h2 to return class 'h2_comparison'.

    # Actually, simpler: just create plot_h2 function? The user asked for "visualize it".
    # Base R 'plot' S3 on the result dataframe is best.

    stop("Use plot.h2_comparison instead")
}
