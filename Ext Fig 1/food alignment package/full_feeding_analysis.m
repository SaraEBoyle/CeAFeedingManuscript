function [] = full_feeding_analysis()
align_to_event_frames();
calculate_z_score_w_consistent_baselinesV2;
remove_high_artifact_trials_V3();
calculate_z_score_w_consistent_baselinesV2;
%plot_events_whole_session;
end