"""Initial DB structure

Revision ID: 2de031b42935
Revises: 
Create Date: 2021-02-09 10:55:01.062321

"""
from alembic import op
import sqlalchemy as sa
import sqlalchemy.dialects.postgresql as pgsql
from sqlalchemy.schema import Sequence, CreateSequence

# revision identifiers, used by Alembic.
revision = '2de031b42935'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    op.execute(CreateSequence(Sequence('molecule_id_seq')))
    op.create_table(
        "molecule",
        sa.Column("molecule_id", sa.Integer(), nullable=False),
        sa.Column("smiles", sa.Text(), nullable=False, unique=True),
        sa.Column("molecular_weight", sa.Float(), nullable=True),
        sa.Column("dft_data", pgsql.JSONB(), nullable=True),
        sa.Column("xtb_data", pgsql.JSONB(), nullable=True),
        sa.Column("ml_data", pgsql.JSONB(), nullable=True)
    )

    op.execute(('alter sequence public.molecule_id_seq '
                'owned by public.molecule.molecule_id;'))
    op.execute(("alter table only public.molecule "
                "alter column molecule_id set default "
                "nextval('public.molecule_id_seq'::regclass);"))
    op.execute(('alter table only public.molecule '
                'add constraint molecule_pkey primary key (molecule_id);'))

    op.execute(CreateSequence(Sequence('conformer_id_seq')))
    op.create_table(
        "conformer",
        sa.Column("conformer_id", sa.Integer(), nullable=False),
        sa.Column("molecule_id", sa.Integer(), sa.ForeignKey('molecule.molecule_id')),
        sa.Column("coords", pgsql.ARRAY(sa.Float()), nullable=False),
        sa.Column("elements", pgsql.ARRAY(sa.String), nullable=False),
        sa.Column("data", pgsql.JSONB(), nullable=True)
    )

    op.execute(('alter sequence public.conformer_id_seq '
                'owned by public.conformer.conformer_id;'))
    op.execute(("alter table only public.conformer "
                "alter column conformer_id set default "
                "nextval('public.conformer_id_seq'::regclass);"))
    op.execute(('alter table only public.conformer '
                'add constraint conformer_pkey '
                'primary key (conformer_id);'))


def downgrade():
    op.drop_table("conformer")
    op.drop_table("molecule")
