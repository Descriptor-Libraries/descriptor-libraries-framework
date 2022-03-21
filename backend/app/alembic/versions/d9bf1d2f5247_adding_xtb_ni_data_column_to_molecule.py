"""Adding xtb_ni_data column to molecule

Revision ID: d9bf1d2f5247
Revises: ab470158c630
Create Date: 2021-05-04 09:12:05.950016

"""
import sqlalchemy as sa
import sqlalchemy.dialects.postgresql as pgsql
from alembic import op

# revision identifiers, used by Alembic.
revision = "d9bf1d2f5247"
down_revision = "ab470158c630"
branch_labels = None
depends_on = None


def upgrade():
    op.add_column("molecule", sa.Column("xtb_ni_data", pgsql.JSONB(), nullable=True))


def downgrade():
    op.drop_column("molecule", "xtb_ni_data")
