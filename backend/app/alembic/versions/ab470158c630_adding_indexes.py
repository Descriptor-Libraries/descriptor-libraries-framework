"""Adding indexes

Revision ID: ab470158c630
Revises: 40c80a5f681c
Create Date: 2021-02-18 15:01:10.833802

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'ab470158c630'
down_revision = '40c80a5f681c'
branch_labels = None
depends_on = None


def upgrade():
    op.execute("create index molecule_mol on public.molecule using gist(mol);")
    op.execute("create index molecule_xtb on public.molecule where xtb_data is not null;")
    op.execute("create index molecule_dft on public.molecule where dft_data is not null;")
    op.execute("create index molecule_ml on public.molecule where ml_data is not null;")


def downgrade():
    op.execute("drop index molecule_mol")
    op.execute("drop index molecule_dft")
    op.execute("drop index molecule_xtb")
    op.execute("drop index molecule_ml")
