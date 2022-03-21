"""Adding rdkit extension

Revision ID: 40c80a5f681c
Revises: 2de031b42935
Create Date: 2021-02-09 11:37:21.339181

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '40c80a5f681c'
down_revision = '2de031b42935'
branch_labels = None
depends_on = None


def upgrade():
    op.execute("create extension if not exists rdkit;")
    op.execute('alter table molecule add column mol mol')
    op.execute("""create or replace function rdk_add_mol()
                  returns trigger as $function$
                  declare
                  begin
                    if NEW.smiles is null then
                        raise exception 'smiles can not be null';
                    end if;
                    if NEW.mol is null then
                        NEW.mol = (select mol_from_smiles(NEW.smiles::cstring));
                    end if;
                    if NEW.molecular_weight is null then
                        NEW.molecular_weight = (select mol_amw(NEW.mol));
                    end if;
                    return NEW;
                  end;
                  $function$
                  language plpgsql;""")

    op.execute(("create trigger rdk_add_mol_trigger "
                "before insert on molecule "
                "for each row execute procedure rdk_add_mol();"))


def downgrade():
    op.execute('drop trigger if exists rdk_add_mol_trigger on molecule')
    op.execute('drop function rdk_add_mol')
    op.execute('alter table molecule drop column mol')
    op.execute('drop extension if exists rdkit')
