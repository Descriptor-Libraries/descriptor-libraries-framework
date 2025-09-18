--
-- PostgreSQL database dump
--

-- Dumped from database version 16.2 (Debian 16.2-1.pgdg120+2)
-- Dumped by pg_dump version 16.2 (Debian 16.2-1.pgdg120+2)

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

--
-- Name: cube; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS cube WITH SCHEMA public;


--
-- Name: EXTENSION cube; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION cube IS 'data type for multidimensional cubes';


--
-- Name: rdkit; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS rdkit WITH SCHEMA public;


--
-- Name: EXTENSION rdkit; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION rdkit IS 'Cheminformatics functionality for PostgreSQL.';


SET default_tablespace = '';

SET default_table_access_method = heap;

--
-- Name: conformer; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.conformer (
    molecule_id text,
    coords double precision[],
    elements character varying[],
    basis text,
    method text,
    program text,
    cluster boolean,
    conformer_id integer NOT NULL
);


ALTER TABLE public.conformer OWNER TO postgres;

--
-- Name: conformer_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.conformer_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE public.conformer_id_seq OWNER TO postgres;

--
-- Name: conformer_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.conformer_id_seq OWNED BY public.conformer.conformer_id;


--
-- Name: dft_data; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.dft_data (
    molecule_id text,
    property text,
    boltz numeric,
    low_e numeric,
    max numeric,
    min numeric
);


ALTER TABLE public.dft_data OWNER TO postgres;

--
-- Name: molecule; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.molecule (
    molecule_id text NOT NULL,
    smiles text,
    canonical_smiles text,
    mol public.mol,
    pca public.cube,
    umap public.cube,
    molecular_weight double precision,
    pat text,
    morganbv public.bfp
);


ALTER TABLE public.molecule OWNER TO postgres;

--
-- Name: molecule_molecule_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.molecule_molecule_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE public.molecule_molecule_id_seq OWNER TO postgres;

--
-- Name: molecule_molecule_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.molecule_molecule_id_seq OWNED BY public.molecule.molecule_id;


--
-- Name: qikprop_data; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.qikprop_data (
    molecule_id text NOT NULL,
    property text NOT NULL,
    value numeric
);


ALTER TABLE public.qikprop_data OWNER TO postgres;

--
-- Name: conformer conformer_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.conformer ALTER COLUMN conformer_id SET DEFAULT nextval('public.conformer_id_seq'::regclass);


--
-- Name: molecule molecule_id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.molecule ALTER COLUMN molecule_id SET DEFAULT nextval('public.molecule_molecule_id_seq'::regclass);


--
-- Name: conformer conformer_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.conformer
    ADD CONSTRAINT conformer_pkey PRIMARY KEY (conformer_id);


--
-- Name: molecule molecule_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.molecule
    ADD CONSTRAINT molecule_pkey PRIMARY KEY (molecule_id);


--
-- Name: molecule molecule_smiles_key; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.molecule
    ADD CONSTRAINT molecule_smiles_key UNIQUE (smiles);


--
-- Name: idx_molecule_molecule_id; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX idx_molecule_molecule_id ON public.molecule USING btree (molecule_id);


--
-- Name: idx_molecule_pca; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX idx_molecule_pca ON public.molecule USING btree (pca);


--
-- Name: idx_molecule_pca_knn; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX idx_molecule_pca_knn ON public.molecule USING gist (pca);


--
-- Name: idx_molecule_smiles; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX idx_molecule_smiles ON public.molecule USING btree (smiles);


--
-- Name: idx_molecule_umap; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX idx_molecule_umap ON public.molecule USING btree (umap);


--
-- Name: idx_molecule_umap_knn; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX idx_molecule_umap_knn ON public.molecule USING gist (umap);


--
-- Name: molecule_mol; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX molecule_mol ON public.molecule USING gist (mol);


--
-- Name: molecule_mol_idx; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX molecule_mol_idx ON public.molecule USING gist (mol);


--
-- Name: morganbv_idx; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX morganbv_idx ON public.molecule USING gist (morganbv);


--
-- Name: mw_idx; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX mw_idx ON public.molecule USING btree (molecular_weight);


--
-- Name: conformer conformer_molecule_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.conformer
    ADD CONSTRAINT conformer_molecule_id_fkey FOREIGN KEY (molecule_id) REFERENCES public.molecule(molecule_id);


--
-- Name: dft_data dft_data_molecule_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.dft_data
    ADD CONSTRAINT dft_data_molecule_id_fkey FOREIGN KEY (molecule_id) REFERENCES public.molecule(molecule_id);


--
-- Name: qikprop_data qikprop_data_molecule_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.qikprop_data
    ADD CONSTRAINT qikprop_data_molecule_id_fkey FOREIGN KEY (molecule_id) REFERENCES public.molecule(molecule_id) ON DELETE CASCADE;


--
-- PostgreSQL database dump complete
--

