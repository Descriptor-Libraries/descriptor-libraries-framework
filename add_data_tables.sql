CREATE TABLE ml_data (
    molecule_id INTEGER,
    property TEXT,
    max DOUBLE PRECISION,
    min DOUBLE PRECISION,
    delta DOUBLE PRECISION,
    vburminconf DOUBLE PRECISION,
    boltzmann_average DOUBLE PRECISION
);

ALTER TABLE ml_data
ADD CONSTRAINT fk_molecule_id
FOREIGN KEY (molecule_id) REFERENCES molecule(molecule_id);

CREATE INDEX idx_ml_data_molecule_id ON ml_data(molecule_id);

\COPY ml_data FROM 'ml_data_json_table.csv' DELIMITER ',' CSV HEADER;

CREATE TABLE dft_data (
    molecule_id INTEGER,
    property TEXT,
    max DOUBLE PRECISION,
    min DOUBLE PRECISION,
    delta DOUBLE PRECISION,
    vburminconf DOUBLE PRECISION,
    boltzmann_average DOUBLE PRECISION
);

ALTER TABLE dft_data
ADD CONSTRAINT fk_molecule_id
FOREIGN KEY (molecule_id) REFERENCES molecule(molecule_id);

CREATE INDEX idx_dft_data_molecule_id ON dft_data(molecule_id);

\COPY dft_data FROM 'dft_data_json_table.csv' DELIMITER ',' CSV HEADER;

CREATE TABLE xtb_data (
    molecule_id INTEGER,
    property TEXT,
    max DOUBLE PRECISION,
    min DOUBLE PRECISION,
    boltzmann_average DOUBLE PRECISION
);

ALTER TABLE xtb_data
ADD CONSTRAINT fk_molecule_id
FOREIGN KEY (molecule_id) REFERENCES molecule(molecule_id);

CREATE INDEX idx_xtb_data_molecule_id ON xtb_data(molecule_id);

\COPY xtb_data FROM 'xtb_data_json_table.csv' DELIMITER ',' CSV HEADER;


CREATE TABLE xtb_ni_data (
    molecule_id INTEGER,
    property TEXT,
    boltzmann_average DOUBLE PRECISION,
    max DOUBLE PRECISION,
    min DOUBLE PRECISION
);

ALTER TABLE xtb_ni_data
ADD CONSTRAINT fk_molecule_id
FOREIGN KEY (molecule_id) REFERENCES molecule(molecule_id);

CREATE INDEX idx_xtb_ni_data_molecule_id ON xtb_ni_data(molecule_id);

\COPY xtb_ni_data FROM 'xtb_ni_data_json_table.csv' DELIMITER ',' CSV HEADER;
