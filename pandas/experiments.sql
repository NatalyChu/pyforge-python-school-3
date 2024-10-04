WITH well_properties AS (
    SELECT well_id, well_row, well_column, plate_id, property_name, property_value
    FROM wells
),
plate_properties AS (
    SELECT plate_id, property_name, property_value
    FROM plates
),
experiment_properties AS (
    SELECT experiment_id, property_name, property_value
    FROM experiments
),
well_plate_properties AS (
    SELECT w.well_id, w.well_row, w.well_column, w.plate_id, 
           COALESCE(w.property_value, p.property_value) AS property_value, 
           w.property_name
    FROM well_properties w
    LEFT JOIN plate_properties p
    ON w.plate_id = p.plate_id AND w.property_name = p.property_name
),
all_properties AS (
    SELECT wpp.well_id, wpp.well_row, wpp.well_column, 
           COALESCE(wpp.property_value, e.property_value) AS property_value, 
           wpp.property_name
    FROM well_plate_properties wpp
    LEFT JOIN experiment_properties e
    ON wpp.experiment_id = e.experiment_id AND wpp.property_name = e.property_name
)
SELECT well_id, well_row, well_column, property_name, property_value
FROM all_properties;
