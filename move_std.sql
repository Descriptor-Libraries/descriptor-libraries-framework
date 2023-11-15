ALTER TABLE xtb_data
ADD COLUMN std DOUBLE PRECISION;

UPDATE xtb_data
SET std = x.boltzmann_average
FROM (
    SELECT 
        SUBSTRING(property FROM 1 FOR POSITION('_std' IN property) - 1) AS property_name, 
        boltzmann_average
    FROM xtb_data
    WHERE property LIKE '%_std%'
) AS x
WHERE xtb_data.property = x.property_name;

DELETE FROM xtb_data
WHERE property LIKE '%_std%';

---------------------------------------------
ALTER TABLE xtb_ni_data
ADD COLUMN std DOUBLE PRECISION;

UPDATE xtb_ni_data
SET std = x.boltzmann_average
FROM (
    SELECT 
        SUBSTRING(property FROM 1 FOR POSITION('_std' IN property) - 1) AS property_name, 
        boltzmann_average
    FROM xtb_ni_data
    WHERE property LIKE '%_std%'
) AS x
WHERE xtb_ni_data.property = x.property_name;

DELETE FROM xtb_ni_data
WHERE property LIKE '%_std%';
