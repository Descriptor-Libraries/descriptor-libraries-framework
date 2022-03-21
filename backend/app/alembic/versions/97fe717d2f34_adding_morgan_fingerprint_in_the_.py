"""adding morgan fingerprint in the database

Revision ID: 97fe717d2f34
Revises: d9bf1d2f5247
Create Date: 2021-05-04 11:10:35.514162

"""
import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision = "97fe717d2f34"
down_revision = "d9bf1d2f5247"
branch_labels = None
depends_on = None


def upgrade():
    op.execute("alter table molecule add column fingerprint bfp")
    op.execute("update molecule set fingerprint = rdkit_fp(mol);")


def downgrade():
    op.execute("alter table molecule drop column fingerprint")
