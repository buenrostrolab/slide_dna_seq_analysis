function[fit_all] = apply_smooth(fit_auto,track_auto,track_all)

    knn = knnsearch(track_auto,track_all,'K',1);
    fit_all = fit_auto(knn);

end